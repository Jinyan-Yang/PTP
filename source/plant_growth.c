/* ============================================================================
* Calls photosynthesis model, water balance and evolves aboveground plant
* C & Nstate. Pools recieve C through allocation of accumulated photosynthate
* and N from both soil uptake and retranslocation within the plant. Key feedback
* through soil N mineralisation and plant N uptake
*
*
* NOTES:
*
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   17.02.2015
*
* =========================================================================== */
#include "plant_growth.h"
#include "water_balance.h"
#include "zbrent.h"
#include <math.h>



void calc_day_growth(canopy_wk *cw, control *c, fluxes *f, fast_spinup *fs,
                     met_arrays *ma, met *m, nrutil *nr, params *p, state *s,
                     double day_length, int doy, double fdecay, double rdecay)
{
    double previous_topsoil_store, dummy=0.0,
           previous_rootzone_store, nitfac, ncbnew, nccnew, ncwimm, ncwnew;
    double previous_sw, current_sw, previous_cs, current_cs, year;
    int    recalc_wb;

    /* Store the previous days soil water store */
    previous_topsoil_store = s->pawater_topsoil;
    previous_rootzone_store = s->pawater_root;

    previous_sw = s->pawater_topsoil + s->pawater_root;
    previous_cs = s->canopy_store;
    year = ma->year[c->day_idx];

    if (c->sub_daily) {
        /* calculate 30 min two-leaf GPP/NPP, respiration and water fluxes */
        canopy(cw, c, f, ma, m, nr, p, s);
    } else {
        /* calculate daily GPP/NPP, respiration and update water balance */
        carbon_daily_production(c, f, m, p, s, day_length);
        calculate_water_balance(c, f, m, p, s, day_length, dummy, dummy, dummy);

        current_sw = s->pawater_topsoil + s->pawater_root;
        current_cs = s->canopy_store;
        f->day_ppt = m->rain;
        //check_water_balance(c, f, s, previous_sw, current_sw, previous_cs,
        //                    current_cs, year, doy);
    }

    /* leaf N:C as a fraction of Ncmaxyoung, i.e. the max N:C ratio of
       foliage in young stand */
    nitfac = MIN(1.0, s->shootnc / p->ncmaxfyoung);

    ///* figure out the C allocation fractions */
    //if (c->deciduous_model){
    //    /* Allocation is annually for deciduous "tree" model, but we need to
    //       keep a check on stresses during the growing season and the LAI
    //       figure out limitations during leaf growth period. This also
    //       applies for deciduous grasses, need to do the growth stress
    //       calc for grasses here too. */
    //    if (s->leaf_out_days[doy] > 0.0) {

    //        calc_carbon_allocation_fracs(c, f, fs, p, s, ma,nitfac);

    //        /* store the days allocation fraction, we average these at the
    //           end of the year (for the growing season) */
    //        s->avg_alleaf += f->alleaf;
    //        s->avg_albranch += f->albranch;
    //        s->avg_alstem += f->alstem;
    //        s->avg_alroot += f->alroot;
    //        s->avg_alcroot += f->alcroot;

    //    }
    //} else {
    //    /* daily allocation...*/
    //    calc_carbon_allocation_fracs(c, f, fs, p, s,ma, nitfac);
    //}
	calc_carbon_allocation_fracs(c, f, fs, p, s, ma, nitfac);
    /* Distribute new C and N through the system */
    carbon_allocation(c, f, p, s, nitfac, doy);

    calculate_ncwood_ratios(c, p, s, nitfac, &ncbnew, &nccnew, &ncwimm,
                            &ncwnew);

    recalc_wb = nitrogen_allocation(c, f, p, s, ncbnew, nccnew, ncwimm, ncwnew,
                                    fdecay, rdecay, doy);

    //if (c->exudation && c->alloc_model != GRASSES) {
    //    calc_root_exudation(c, f, p, s);
    //}

    /* If we didn't have enough N available to satisfy wood demand, NPP
       is down-regulated and thus so is GPP. We also need to recalculate the
       water balance given the lower GPP. */
    if (recalc_wb) {
        s->pawater_topsoil = previous_topsoil_store;
        s->pawater_root = previous_rootzone_store;

        if (c->sub_daily) {
            /* reduce transpiration to match cut back GPP
                -there isn't an obvious way to make this work at the 30 min
                 timestep, so invert T from WUE assumption and use that
                 to recalculate the end day water balance
            */
            f->transpiration = f->gpp_gCm2 / f->wue;
            update_water_storage_recalwb(c, f, p, s, m);

        } else {
            calculate_water_balance(c, f, m, p, s, day_length, dummy, dummy,
                                    dummy);
        }

    }
    update_plant_state(c, f, p, s, fdecay, rdecay, doy);

    precision_control(f, s);

    return;
}

//allocation fraction
void calc_carbon_allocation_fracs(control* c, fluxes* f, fast_spinup* fs,
	params* p, state* s, met_arrays* ma, double nitfac) {
	/* Carbon allocation fractions to move photosynthate through the plant.

	Parameters:
	-----------
	nitfac : float
		leaf N:C as a fraction of 'Ncmaxfyoung' (max 1.0)

	Returns:
	--------
	alleaf : float
		allocation fraction for shoot
	alroot : float
		allocation fraction for fine roots
	albranch : float
		allocation fraction for branches
	alstem : float
		allocation fraction for stem

	References:
	-----------
	Corbeels, M. et al (2005) Ecological Modelling, 187, 449-474.
	McMurtrie, R. E. et al (2000) Plant and Soil, 224, 135-152.
	*/
	double min_leaf_alloc, adj, arg1, arg2, arg3, arg4, leaf2sa_target,
		sap_cross_sec_area, lr_max, stress, mis_match, orig_ar,
		reduction, target_branch, coarse_root_target, left_over,
		total_alloc, leaf2sap, spare, gLeaf, gRoot, ppt_sum_prev,
		t_et, plant_cover;
	int dd, pass_day;
	/* this is obviously arbitary */
	double min_stem_alloc = 0.01;

	if (c->alloc_model == FIXED) {
		f->alleaf = (p->c_alloc_fmax + nitfac *
			(p->c_alloc_fmax - p->c_alloc_fmin));

		f->alroot = (p->c_alloc_rmax + nitfac *
			(p->c_alloc_rmax - p->c_alloc_rmin));

		f->albranch = (p->c_alloc_bmax + nitfac *
			(p->c_alloc_bmax - p->c_alloc_bmin));

		/* allocate remainder to stem */
		f->alstem = 1.0 - f->alleaf - f->alroot - f->albranch;

		f->alcroot = p->c_alloc_cmax * f->alstem;
		f->alstem -= f->alcroot;

	}
    else if (c->alloc_model == GRASSES) {
        /* First figure out root allocation given available water & nutrients
              hyperbola shape to allocation */
        f->alroot = (p->c_alloc_rmax * p->c_alloc_rmin /
            (p->c_alloc_rmin + (p->c_alloc_rmax - p->c_alloc_rmin) *
                s->prev_sma));
        f->alleaf = 1.0 - f->alroot;

        /* Now adjust root & leaf allocation to maintain balance, accounting
           for stress e.g. -> Sitch et al. 2003, GCB.
         leaf-to-root ratio under non-stressed conditons
        lr_max = 0.8;
         Calculate adjustment on lr_max, based on current "stress"
           calculated from running mean of N and water stress
        stress = lr_max * s->prev_sma;
        calculate new allocation fractions based on imbalance in *biomass*
        mis_match = s->shoot / (s->root * stress);
        if (mis_match > 1.0) {
            reduce leaf allocation fraction
            adj = f->alleaf / mis_match;
            f->alleaf = MAX(p->c_alloc_fmin, MIN(p->c_alloc_fmax, adj));
            f->alroot = 1.0 - f->alleaf;
        } else {
             reduce root allocation
            adj = f->alroot * mis_match;
            f->alroot = MAX(p->c_alloc_rmin, MIN(p->c_alloc_rmax, adj));
            f->alleaf = 1.0 - f->alroot;
        }*/
        f->alstem = 0.0;
        f->albranch = 0.0;
        f->alcroot = 0.0;

    }
	else if (c->alloc_model == ALLOMETRIC) {

		/* Calculate tree height: allometric reln using the power function
		   (Causton, 1985) */
		s->canht = p->heighto * pow(s->stem, p->htpower);

		/* LAI to stem sapwood cross-sectional area (As m-2 m-2)
		   (dimensionless)
		   Assume it varies between LS0 and LS1 as a linear function of tree
		   height (m) */
		arg1 = s->sapwood * TONNES_AS_KG * M2_AS_HA;
		arg2 = s->canht * p->density * p->cfracts;
		sap_cross_sec_area = arg1 / arg2;
		leaf2sap = s->lai / sap_cross_sec_area;

		/* Allocation to leaves dependant on height. Modification of pipe
		   theory, leaf-to-sapwood ratio is not constant above a certain
		   height, due to hydraulic constraints (Magnani et al 2000; Deckmyn
		   et al. 2006). */

		if (s->canht < p->height0) {
			leaf2sa_target = p->leafsap0;
		}
		else if (float_eq(s->canht, p->height1)) {
			leaf2sa_target = p->leafsap1;
		}
		else if (s->canht > p->height1) {
			leaf2sa_target = p->leafsap1;
		}
		else {
			arg1 = p->leafsap0;
			arg2 = p->leafsap1 - p->leafsap0;
			arg3 = s->canht - p->height0;
			arg4 = p->height1 - p->height0;
			leaf2sa_target = arg1 + (arg2 * arg3 / arg4);
		}
		f->alleaf = alloc_goal_seek(leaf2sap, leaf2sa_target, p->c_alloc_fmax,
			p->targ_sens);

		/* Allocation to branch dependent on relationship between the stem
		   and branch */
		target_branch = p->branch0 * pow(s->stem, p->branch1);
		f->albranch = alloc_goal_seek(s->branch, target_branch, p->c_alloc_bmax,
			p->targ_sens);

		coarse_root_target = p->croot0 * pow(s->stem, p->croot1);
		f->alcroot = alloc_goal_seek(s->croot, coarse_root_target,
			p->c_alloc_cmax, p->targ_sens);

		/* figure out root allocation given available water & nutrients
		   hyperbola shape to allocation, this is adjusted below as we aim
		   to maintain a functional balance */

		f->alroot = (p->c_alloc_rmax * p->c_alloc_rmin /
			(p->c_alloc_rmin + (p->c_alloc_rmax - p->c_alloc_rmin) *
				s->prev_sma));

		f->alstem = 1.0 - f->alroot - f->albranch - f->alleaf - f->alcroot;

		/* minimum allocation to leaves - without it tree would die, as this
		   is done annually. */
		if (c->deciduous_model) {
			if (f->alleaf < 0.05) {
				min_leaf_alloc = 0.05;
				if (f->alstem > min_leaf_alloc)
					f->alstem -= min_leaf_alloc;
				else
					f->alroot -= min_leaf_alloc;
				f->alleaf = min_leaf_alloc;
			}
		}
	}
	else if (c->alloc_model == SGS) {
		//the SGS scheme
		t_et = s->wtfac_topsoil;//f->transpiration / (f->et * s->fipar);

		if (t_et > 1) {
			t_et = 1.0;
		}
		//printf("t-et: %f\n", t_et);
		f->alleaf = t_et * p->c_alloc_fmax;
        //this assume that plant can store c!
        //however one need to make sure that froot and fleaf sum is less than 1
		//f->alroot = (p->c_alloc_rmax + p->c_alloc_fmax) - f->alleaf;
        f->alroot = p->c_alloc_rmax;
	}
	else if(c->alloc_model == FATICHI){
	////calculate past rainfall
	//ppt_sum_prev = 0.0;
	//pass_day = 0;
	//for (dd = c->day_idx;
	//	(dd - pass_day) > 1 && pass_day < p->days_rain;
	//	pass_day++) {

	//	if (pass_day > 0) {
	//		ppt_sum_prev += ma->rain[dd - pass_day];
	//	}
	//	else {
	//		ppt_sum_prev = 0.0;
	//	}
	//}
	//
	//	//	a rough approximation of Hufkens model	
	//if (ppt_sum_prev > p->green_sw_frac) {
	//	gLeaf = p->c_alloc_fmax;
	//	gRoot = p->c_alloc_rmax;
	//}
	//else {
	//	gLeaf = 0.0;
	//	gRoot = 0.0;
	//}

	//f->alroot = gRoot * (1 - s->wtfac_root);
	//f->alleaf = gLeaf * s->wtfac_root;

	//ppt_sum_prev = 0.0;
	//pass_day = 0;

	//f->alstem = 0.0;
	//f->albranch = 0.0;
	//f->alcroot = 0.0;
}
    else if (c->alloc_model == HUFKEN) {

    ///* First figure out root allocation given available water & nutrients
    //   hyperbola shape to allocation */
    //f->alroot = (p->c_alloc_rmax * p->c_alloc_rmin /
    //             (p->c_alloc_rmin + (p->c_alloc_rmax - p->c_alloc_rmin) *
    //              s->prev_sma));
    //f->alleaf = 1.0 - f->alroot;

    //add new allocation based on water
    //calculate past rainfall

    ppt_sum_prev = 0.0;
    pass_day = 0;
    for (dd = c->day_idx;
        (dd - pass_day) > 1 && pass_day < p->days_rain;
        pass_day++) {

        if (pass_day > 0) {
            ppt_sum_prev += ma->rain[dd - pass_day];
        }
        else {
            ppt_sum_prev = 0.0;
        }
    }

    /*
        if (s->wtfac_root > p->green_sw_frac) {
            gLeaf = p->c_alloc_fmax;
            gRoot = p->c_alloc_rmax;
        }
        else {
            gLeaf = 0.0;
            gRoot = 0.0;
        }
*/
//	a rough approximation of Hufkens model	
    if (ppt_sum_prev > p->green_sw_frac) {
        gLeaf = p->c_alloc_fmax;
        gRoot = p->c_alloc_rmax;
    }
    else {
        gLeaf = 0.0;
        gRoot = 0.0;
    }

    // gLeaf = p->c_alloc_fmax;
     //gRoot = p->c_alloc_rmax;
 

    //f->alroot = gRoot * (1 - s->wtfac_topsoil);
    //f->alroot = p->c_alloc_rmax;
    // reduce growth when close to full cover; 
 //here I used the apar function; 
 //maybe better witha a s shape one
    plant_cover = (1 - exp(-0.5 * s->lai))*p->use_cover;
    f->alleaf = gLeaf * pow(s->wtfac_topsoil, p->q) * (1 - plant_cover);
    f->alroot = 1 - f->alleaf;
    ppt_sum_prev = 0.0;
    pass_day = 0;
    /* Now adjust root & leaf allocation to maintain balance, accounting
       for stress e.g. -> Sitch et al. 2003, GCB.

     leaf-to-root ratio under non-stressed conditons
    lr_max = 0.8;

     Calculate adjustment on lr_max, based on current "stress"
       calculated from running mean of N and water stress
    stress = lr_max * s->prev_sma;

    calculate new allocation fractions based on imbalance in *biomass*
    mis_match = s->shoot / (s->root * stress);


    if (mis_match > 1.0) {
        reduce leaf allocation fraction
        adj = f->alleaf / mis_match;
        f->alleaf = MAX(p->c_alloc_fmin, MIN(p->c_alloc_fmax, adj));
        f->alroot = 1.0 - f->alleaf;
    } else {
         reduce root allocation
        adj = f->alroot * mis_match;
        f->alroot = MAX(p->c_alloc_rmin, MIN(p->c_alloc_rmax, adj));
        f->alleaf = 1.0 - f->alroot;
    }*/
    f->alstem = 0.0;
    f->albranch = 0.0;
    f->alcroot = 0.0;

    }
	else {
		fprintf(stderr, "Unknown C allocation model: %d\n", c->alloc_model);
		exit(EXIT_FAILURE);
	}

	///*printf("%f %f %f %f %f\n", f->alleaf, f->albranch + f->alstem, f->alroot,  f->alcroot, s->canht);*/

	///* Total allocation should be one, if not print warning */
	total_alloc = f->alroot + f->alleaf + f->albranch + f->alstem + f->alcroot;
	if (total_alloc > 1.0 + EPSILON) {
		fprintf(stderr, "Allocation fracs > 1: %.13f\n", total_alloc);
		exit(EXIT_FAILURE);
	}

	//if (c->spinup_method == SAS) {
	//	fs->alloc[AF] += f->alleaf;
	//	fs->alloc[AR] += f->alroot;
	//	fs->alloc[ACR] += f->alcroot;
	//	fs->alloc[AB] += f->albranch;
	//	fs->alloc[AW] += f->alstem;
	//}

	return;
}
// allocation of carbon
void carbon_allocation(control* c, fluxes* f, params* p, state* s,
	double nitfac, int doy) {
	/* C distribution - allocate available C through system

	Parameters:
	-----------
	nitfac : float
		leaf N:C as a fraction of 'Ncmaxfyoung' (max 1.0)
	*/
	s->nsc += f->npp;
	f->cpleaf = s->nsc * f->alleaf;
	f->cproot = s->nsc * f->alroot;

	s->nsc += -(f->cpleaf + f->cproot);
	s->nsc *= 0.99;// assuming that 1% will lost via respiration

	//double days_left;
	//if (c->deciduous_model) {
	//	days_left = s->growing_days[doy];
	//	f->cpleaf = f->lrate * days_left;
	//	f->cpbranch = f->brate * days_left;
	//	f->cpstem = f->wrate * days_left;
	//	f->cproot = s->c_to_alloc_root * 1.0 / c->num_days;
	//	f->cpcroot = f->crate * days_left;
	//}
	//else {
	//	s->nsc += f->npp;
	//	f->cpleaf = s->nsc * f->alleaf;
	//	f->cproot = s->nsc * f->alroot;

	//	s->nsc -= (f->cpleaf + f->cproot);
	//	s->nsc *= 0.99;// assuming that 1% will lost via respiration

	//	//// this assumes a minimum nsc to ensure plant reemerge after drought
	//	//if (s->nsc < 1.0) {
	//	//	s->nsc += (f->cpleaf + f->cproot);
	//	//	f->cpleaf = 0.0;
	//	//	f->cproot = 0.0;
	//	//}
	//}



	/* evaluate SLA of new foliage accounting for variation in SLA
	   with tree and leaf age (Sands and Landsberg, 2002). Assume
	   SLA of new foliage is linearly related to leaf N:C ratio
	   via nitfac. Based on date from two E.globulus stands in SW Aus, see
	   Corbeels et al (2005) Ecological Modelling, 187, 449-474.
	   (m2 onesided/kg DW) */
	   //p->sla = p->slazero + nitfac * (p->slamax - p->slazero);

	  // if (c->deciduous_model) {
	  //     if (float_eq(s->shoot, 0.0)) {
	  //         s->lai = 0.0;
	  //     } else if (s->leaf_out_days[doy] > 0.0) {
	  //         s->lai += (f->cpleaf *
	  //                   (p->sla * M2_AS_HA / (KG_AS_TONNES * p->cfracts)) -
	  //                   (f->deadleaves + f->ceaten) * s->lai / s->shoot);
	  //     } else {
	  //         s->lai = 0.0;
	  //     }
	  // } else {
	  //     /* update leaf area [m2 m-2] */
	  //     if (float_eq(s->shoot, 0.0)) {
	  //         s->lai = 0.0;
	  //     } else {
	  //         //s->lai += (f->cpleaf *
	  //         //          (p->sla * M2_AS_HA / (KG_AS_TONNES * p->cfracts)) -
	  //         //          (f->deadleaves + f->ceaten) * s->lai / s->shoot);

			   ////s->lai += (f->cpleaf - f->deadleaves - f->ceaten) *
			   ////	p->sla * M2_AS_HA / (KG_AS_TONNES * p->cfracts);
			   //s->lai += s->shoot * (p->sla * M2_AS_HA / (KG_AS_TONNES * p->cfracts));
	  //     }
	  // }

	   // moved the calculation of lai to up date plant state function
	   // this is to ensure  lai is based on the correct leaf biomass

	return;
}
// change in c pool
void update_plant_state(control* c, fluxes* f, params* p, state* s,
	double fdecay, double rdecay, int doy) {
	/*
	Daily change in C content

	Parameters:
	-----------
	fdecay : float
		foliage decay rate
	rdecay : float
		fine root decay rate

	*/

	double age_effect, ncmaxf, ncmaxr, extras, extrar,f_decay_actual;

	/*
	** Carbon pools
	*/

	/*
     if (f->cpleaf > 0.0) {
		s->shoot += f->cpleaf - f->ceaten;
	}
	else {
		s->shoot += -f->deadleaves - f->ceaten;
	}
    */

    //account for now q_s conditions
    /*
    /double soil_water_impact;
    double soil_water_deficit;

        if (float_eq(p->q_s, 0.0)) {
            soil_water_impact = 1.0;
        }
        else {
            soil_water_deficit = (1 - s->wtfac_topsoil);
            soil_water_impact = pow(soil_water_deficit, p->q_s);
        }
        */

    f_decay_actual = p->fdecay * pow((1.0 - s->wtfac_root), p->q_s) * (1.0 - exp(-0.5 * s->lai));

    if (f_decay_actual < 0.002) {
        f_decay_actual = 0.002; 
    }
    
    f->deadleaves = f_decay_actual * s->shoot;

    s->shoot += f->cpleaf - f->deadleaves - f->ceaten;

    if (s->shoot < 0.0 ) {
        s->shoot = 0.0;
    }
    //s->shoot -= p->fdecay * (1 - s->wtfac_topsoil) * s->shoot;
	//s->shoot *= 0.5;
	s->root += f->cproot - f->deadroots;
	//s->root *= 0.5;

	if (float_eq(s->shoot, 0.0)) {
		s->lai = 0.0;
	}
	else {
		s->lai = s->shoot * (p->sla * M2_AS_HA / (KG_AS_TONNES * p->cfracts));
	}

	if (c->fixed_lai) {
		s->lai = p->fix_lai;
	}

	s->croot += f->cpcroot - f->deadcroots;
	s->branch += f->cpbranch - f->deadbranch;
	s->stem += f->cpstem - f->deadstems;

	/* annoying but can't see an easier way with the code as it is.
	   If we are modelling grases, i.e. no stem them without this
	   the sapwood will end up being reduced to a silly number as
	   deadsapwood will keep being removed from the pool, even though there
	   is no wood. */
	if (float_eq(s->stem, 0.01)) {
		s->sapwood = 0.01;
	}
	else if (s->stem < 0.01) {
		s->sapwood = 0.01;
	}
	else {
		s->sapwood += f->cpstem - f->deadsapwood;
	}


	/*
	** Nitrogen pools
	*/
	if (c->deciduous_model) {
		s->shootn += (f->npleaf - (f->lnrate * s->remaining_days[doy]) -
			f->neaten);
	}
	else {
		s->shootn += f->npleaf - fdecay * s->shootn - f->neaten;
	}

	s->branchn += f->npbranch - p->bdecay * s->branchn;
	s->rootn += f->nproot - rdecay * s->rootn;
	s->crootn += f->npcroot - p->crdecay * s->crootn;
	s->stemnimm += f->npstemimm - p->wdecay * s->stemnimm;
	s->stemnmob += (f->npstemmob - p->wdecay * s->stemnmob - p->retransmob *
		s->stemnmob);
	s->stemn = s->stemnimm + s->stemnmob;


	if (c->deciduous_model == FALSE) {
		/*
		   =============================
			Enforce maximum N:C ratios.
		   =============================
		*/

		/* If foliage or root N/C exceeds its max, then N uptake is cut back*/

		/* maximum leaf n:c ratio is function of stand age
			- switch off age effect by setting ncmaxfyoung = ncmaxfold */
		age_effect = (s->age - p->ageyoung) / (p->ageold - p->ageyoung);
		ncmaxf = p->ncmaxfyoung - (p->ncmaxfyoung - p->ncmaxfold) * age_effect;

		if (ncmaxf < p->ncmaxfold)
			ncmaxf = p->ncmaxfold;

		if (ncmaxf > p->ncmaxfyoung)
			ncmaxf = p->ncmaxfyoung;

		// if shoot N:C ratio exceeds its max, then nitrogen uptake is cut back
		extras = 0.0;
		if (s->lai > 0.0) {
			if (s->shootn > (s->shoot * ncmaxf)) {
				extras = s->shootn - s->shoot * ncmaxf;

				/* Ensure N uptake cannot be reduced below zero. */
				if (extras > f->nuptake) {
					extras = f->nuptake;
				}

				s->shootn -= extras;
				f->nuptake -= extras;
			}
		}

		// if root N:C ratio exceeds its max, then nitrogen uptake is cut back

		// max root n:c
		ncmaxr = ncmaxf * p->ncrfac;
		extrar = 0.0;
		if (s->rootn > (s->root * ncmaxr)) {
			extrar = s->rootn - s->root * ncmaxr;

			/* Ensure N uptake cannot be reduced below zero. */
			if (extrar > f->nuptake) {
				extrar = f->nuptake;
			}

			s->rootn -= extrar;
			f->nuptake -= extrar;
		}
	}

	/* Update deciduous storage pools */
	if (c->deciduous_model)
		calculate_cn_store(c, f, s);

	return;
}


void calc_root_exudation(control *c, fluxes *f, params *p, state *s) {
    /*
        Rhizodeposition (f->root_exc) is assumed to be a fraction of the
        current root growth rate (f->cproot), which increases with increasing
        N stress of the plant.
    */
    double CN_leaf, frac_to_rexc, CN_ref, arg;

    if (float_eq(s->shoot, 0.0) || float_eq(s->shootn, 0.0)) {
        /* nothing happens during leaf off period */
        CN_leaf = 0.0;
        frac_to_rexc = 0.0;
    } else {

        if (c->deciduous_model) {
            /* broadleaf */
            CN_ref = 25.0;
        } else {
            /* conifer */
            CN_ref = 42.0;
        }

        /*
        ** The fraction of growth allocated to rhizodeposition, constrained
        ** to solutions lower than 0.5
        */
        CN_leaf = 1.0 / s->shootnc;
        arg = MAX(0.0, (CN_leaf - CN_ref) / CN_ref);
        frac_to_rexc = MIN(0.5, p->a0rhizo + p->a1rhizo * arg);
    }

    /* Rhizodeposition */
    f->root_exc = frac_to_rexc * f->cproot;
    if (float_eq(f->cproot, 0.0)) {
        f->root_exn = 0.0;
    } else {
        /*
        ** N flux associated with rhizodeposition is based on the assumption
        ** that the CN ratio of rhizodeposition is equal to that of fine root
        ** growth
        */
        f->root_exn = f->root_exc * (f->nproot / f->cproot);
    }

    /*
    ** Need to remove exudation C & N fluxes from fine root growth fluxes so
    ** that things balance.
    */
    f->cproot -= f->root_exc;
    f->nproot -= f->root_exn;

    return;
}

void carbon_daily_production(control *c, fluxes *f, met *m, params *p, state *s,
                             double daylen) {
    /* Calculate GPP, NPP and plant respiration at the daily timestep

    Parameters:
    -----------
    daylen : float
        daytime length (hrs)

    References:
    -----------
    * Jackson, J. E. and Palmer, J. W. (1981) Annals of Botany, 47, 561-565.
    */
    double leafn, fc, ncontent;

    if (s->lai > 0.0) {
        /* average leaf nitrogen content (g N m-2 leaf) */
        leafn = (s->shootnc * p->cfracts / p->sla * KG_AS_G);

        /* total nitrogen content of the canopy */
        ncontent = leafn * s->lai;

    } else {
        ncontent = 0.0;
    }

    /* When canopy is not closed, canopy light interception is reduced
        - calculate the fractional ground cover */
    if (s->lai < p->lai_closed) {
        /* discontinuous canopies */
        fc = s->lai / p->lai_closed;
    } else {
        fc = 1.0;
    }

    /* fIPAR - the fraction of intercepted PAR = IPAR/PAR incident at the
       top of the canopy, accounting for partial closure based on Jackson
       and Palmer (1979). */
    if (s->lai > 0.0)
        s->fipar = ((1.0 - exp(-p->kext * s->lai / fc)) * fc);
    else
        s->fipar = 0.0;

    if (c->water_stress) {
        /* Calculate the soil moisture availability factors [0,1] in the
           topsoil and the entire root zone */
        calculate_soil_water_fac(c, p, s);
    } else {
        /* really this should only be a debugging option! */
        s->wtfac_topsoil = 1.0;
        s->wtfac_root = 1.0;
    }
    /* Estimate photosynthesis */
    if (c->assim_model == BEWDY){
        exit(EXIT_FAILURE);
    } else if (c->assim_model == MATE) {
        if (c->ps_pathway == C3) {
            mate_C3_photosynthesis(c, f, m, p, s, daylen, ncontent);
        } else {
            mate_C4_photosynthesis(c, f, m, p, s, daylen, ncontent);
        }
    } else {
        fprintf(stderr,"Unknown photosynthesis model'");
        exit(EXIT_FAILURE);
    }

    /* Calculate plant respiration */
    if (c->respiration_model == FIXED) {
        /* Plant respiration assuming carbon-use efficiency. */
        f->auto_resp = f->gpp * p->cue;
    } else if (c->respiration_model == VARY) {
        calc_autotrophic_respiration(c, f, m, p, s);
    }

    f->npp = MAX(0.0, f->gpp - f->auto_resp);
    f->npp_gCm2 = f->npp * TONNES_HA_2_G_M2;

    return;
}

void calc_autotrophic_respiration(control *c, fluxes *f, met *m, params *p,
                                  state *s) {
    // Autotrophic respiration is the sum of the growth component (Rg)
    // and the the temperature-dependent maintenance respiration (Rm) of
    // leaves (Rml), fine roots (Rmr) and wood (Rmw)

    double Rml, Rmw, Rmr, Rm, Rg, rk, k = 0.0548;
    double shootn, rootn, stemn, cue;

    // respiration rate (gC gN-1 d-1) on a 10degC base
    rk = p->resp_coeff * k;

    if (c->ncycle == FALSE) {
        shootn = (s->shoot * 0.03) * TONNES_HA_2_G_M2;
        rootn = (s->root * 0.02) * TONNES_HA_2_G_M2;
        stemn = (s->stem * 0.003) * TONNES_HA_2_G_M2;
    } else {
        shootn = s->shootn * TONNES_HA_2_G_M2;
        rootn = s->rootn * TONNES_HA_2_G_M2;
        stemn = s->stemn * TONNES_HA_2_G_M2;
    }

    // Maintenance respiration, the cost of metabolic processes in
    // living tissues, differs according to tissue N
    Rml = rk * shootn * lloyd_and_taylor(m->tair);
    Rmw = rk * stemn * lloyd_and_taylor(m->tair); // should really be sapwood
    Rmr = rk * rootn * lloyd_and_taylor(m->tsoil);
    Rm = (Rml + Rmw + Rmr) * GRAM_C_2_TONNES_HA;

    // After maintenance respiration is subtracted from GPP, 25% of the
    // remainder is taken as growth respiration, the cost of producing new
    // tissues
    Rg = MAX(0.0, (f->gpp - Rm) * 0.25);
    f->auto_resp = Rm + Rg;

    // Should be revisited, but it occurs to me that during spinup
    // the tissue initialisation could be greater than the incoming gpp and so
    // we might never grow if we respire all out C. Clearly were we using a
    // storage pool this wouldn't be such a drama. For now bound it...
    //
    // De Lucia et al. Global Change Biology (2007) 13, 1157–1167:
    // CUE varied from 0.23 to 0.83 from a literature survey
    cue = (f->gpp - f->auto_resp) / f->gpp;
    if (cue < 0.2) {
        f->auto_resp = f->gpp * 0.8;
    } else if (cue > 0.8) {
        f->auto_resp = f->gpp * 0.2;
    }

    return;
}

double lloyd_and_taylor(double temp) {
    // Modified Arrhenius equation (Lloyd & Taylor, 1994)
    // The modification introduced by Lloyd & Taylor (1994) represents a
    // decline in the parameter for activation energy with temperature.
    //
    // Parameters:
    // -----------
    // temp : float
    //      temp deg C

    return (exp(308.56 * ((1.0 / 56.02) - (1.0 / (temp + 46.02)))));
}

void calculate_ncwood_ratios(control *c, params *p, state *s, double nitfac,
                             double *ncbnew, double *nccnew, double *ncwimm,
                             double *ncwnew) {
    /* Estimate the N:C ratio in the branch and stem. Option to vary
    the N:C ratio of the stem following Jeffreys (1999) or keep it a fixed
    fraction

    Parameters:
    -----------
    nitfac : float
        leaf N:C as a fraction of the max N:C ratio of foliage in young
        stand

    Returns:
    --------
    ncbnew : float
        N:C ratio of branch
    nccnew : double
        N:C ratio of coarse root
    ncwimm : float
        N:C ratio of immobile stem
    ncwnew : float
        N:C ratio of mobile stem

    References:
    ----------
    * Jeffreys, M. P. (1999) Dynamics of stemwood nitrogen in Pinus radiata
      with modelled implications for forest productivity under elevated
      atmospheric carbon dioxide. PhD.
    */

    /* n:c ratio of new branch wood*/
    *ncbnew = p->ncbnew + nitfac * (p->ncbnew - p->ncbnewz);

    /* n:c ratio of coarse root */
    *nccnew = p->nccnew + nitfac * (p->nccnew - p->nccnewz);

    /* fixed N:C in the stemwood */
    if (c->fixed_stem_nc) {
        /* n:c ratio of stemwood - immobile pool and new ring */
        *ncwimm = p->ncwimm + nitfac * (p->ncwimm - p->ncwimmz);

        /* New stem ring N:C at critical leaf N:C (mobile) */
        *ncwnew = p->ncwnew + nitfac * (p->ncwnew - p->ncwnewz);

   /* vary stem N:C based on reln with foliage, see Jeffreys. Jeffreys 1999
      showed that N:C ratio of new wood increases with foliar N:C ratio,
      modelled here based on evidence as a linear function. */
   } else {
        *ncwimm = MAX(0.0, (0.0282 * s->shootnc + 0.000234) * p->fhw);

        /* New stem ring N:C at critical leaf N:C (mobile) */
        *ncwnew = MAX(0.0, 0.162 * s->shootnc - 0.00143);
    }

    return;
}


// nitrogen bit
int nitrogen_allocation(control *c, fluxes *f, params *p, state *s,
                        double ncbnew, double nccnew, double ncwimm,
                        double ncwnew, double fdecay, double rdecay, int doy) {
    /* Nitrogen distribution - allocate available N through system.
    N is first allocated to the woody component, surplus N is then allocated
    to the shoot and roots with flexible ratios.

    References:
    -----------
    McMurtrie, R. E. et al (2000) Plant and Soil, 224, 135-152.

    Parameters:
    -----------
    ncbnew : float
        N:C ratio of branch
    ncwimm : float
        N:C ratio of immobile stem
    ncwnew : float
        N:C ratio of mobile stem
    fdecay : float
        foliage decay rate
    rdecay : float
        fine root decay rate
    */

    int    recalc_wb;
    double nsupply, rtot, ntot, arg, lai_inc = 0.0, conv;
    double depth_guess = 1.0;

    /* default is we don't need to recalculate the water balance,
       however if we cut back on NPP due to available N below then we do
       need to do this */
    recalc_wb = FALSE;

    /* N retranslocated proportion from dying plant tissue and stored within
       the plant */
    f->retrans = nitrogen_retrans(c, f, p, s, fdecay, rdecay, doy);
    f->nuptake = calculate_nuptake(c, p, s);

    /*  Ross's Root Model. */
    if (c->model_optroot) {

        /* convert t ha-1 day-1 to gN m-2 year-1 */
        nsupply = (calculate_nuptake(c, p, s) *
                   TONNES_HA_2_G_M2 * DAYS_IN_YRS);

        /* covnert t ha-1 to kg DM m-2 */
        rtot = s->root * TONNES_HA_2_KG_M2 / p->cfracts;
        /*f->nuptake_old = f->nuptake; */

        calc_opt_root_depth(p->d0x, p->r0, p->topsoil_depth * MM_TO_M,
                            rtot, nsupply, depth_guess, &s->root_depth,
                            &f->nuptake, &f->rabove);

        /*umax = self.rm.calc_umax(f->nuptake) */

        /* covert nuptake from gN m-2 year-1  to t ha-1 day-1 */
        f->nuptake = f->nuptake * G_M2_2_TONNES_HA * YRS_IN_DAYS;

        /* covert from kg DM N m-2 to t ha-1 */
        f->deadroots = p->rdecay * f->rabove * p->cfracts * KG_M2_2_TONNES_HA;
        f->deadrootn = s->rootnc * (1.0 - p->rretrans) * f->deadroots;
    }

    /* Mineralised nitrogen lost from the system by volatilisation/leaching */
    f->nloss = p->rateloss * s->inorgn;

    /* total nitrogen to allocate */
    ntot = MAX(0.0, f->nuptake + f->retrans);

    if (c->deciduous_model) {
        /* allocate N to pools with fixed N:C ratios */

        /* N flux into new ring (immobile component -> structrual components) */
        f->npstemimm = f->wnimrate * s->growing_days[doy];

        /* N flux into new ring (mobile component -> can be retrans for new
           woody tissue) */
        f->npstemmob = f->wnmobrate * s->growing_days[doy];
        f->nproot = s->n_to_alloc_root / c->num_days;
        f->npcroot = f->cnrate * s->growing_days[doy];
        f->npleaf = f->lnrate * s->growing_days[doy];
        f->npbranch = f->bnrate * s->growing_days[doy];
    } else {
        /* allocate N to pools with fixed N:C ratios */

        /* N flux into new ring (immobile component -> structural components) */
        f->npstemimm = f->npp * f->alstem * ncwimm;

        /* N flux into new ring (mobile component -> can be retrans for new
           woody tissue) */
        f->npstemmob = f->npp * f->alstem * (ncwnew - ncwimm);
        f->npbranch = f->npp * f->albranch * ncbnew;
        f->npcroot = f->npp * f->alcroot * nccnew;

        /* If we have allocated more N than we have available
            - cut back C prodn */
        arg = f->npstemimm + f->npstemmob + f->npbranch + f->npcroot;


        if (arg > ntot && c->fixleafnc == FALSE && c->fixed_lai && c->ncycle) {


            /* Need to readjust the LAI for the reduced growth as this will
               have already been increased. First we need to figure out how
               much we have increased LAI by, important it is done here
               before cpleaf is reduced! */
            if (float_eq(s->shoot, 0.0)) {
                lai_inc = 0.0;
            } else {
                lai_inc = (f->cpleaf *
                           (p->sla * M2_AS_HA / (KG_AS_TONNES * p->cfracts)) -
                           (f->deadleaves + f->ceaten) * s->lai / s->shoot);
            }

            f->npp *= ntot / (f->npstemimm + f->npstemmob + \
                              f->npbranch + f->npcroot);

            /* need to adjust growth values accordingly as well */
			s->nsc += f->npp;
			f->cpleaf = s->nsc * f->alleaf;
            f->cproot = s->nsc * f->alroot;
            f->cpcroot = s->nsc * f->alcroot;
            f->cpbranch = s->nsc * f->albranch;
            f->cpstem = s->nsc * f->alstem;
			s->nsc += -(f->cpleaf - f->cproot - f->cpcroot -
				f->cpbranch - f->cpstem);

            f->npbranch = f->npp * f->albranch * ncbnew;
            f->npstemimm = f->npp * f->alstem * ncwimm;
            f->npstemmob = f->npp * f->alstem * (ncwnew - ncwimm);
            f->npcroot = f->npp * f->alcroot * nccnew;

            /* Save WUE before cut back */
            f->wue = f->gpp_gCm2 / f->transpiration;

            /* Also need to recalculate GPP and thus Ra and return a flag
               so that we know to recalculate the water balance. */
            f->gpp = f->npp / p->cue;
            conv = G_AS_TONNES / M2_AS_HA;
            f->gpp_gCm2 = f->gpp / conv;
            f->gpp_am = f->gpp_gCm2 / 2.0;
            f->gpp_pm = f->gpp_gCm2 / 2.0;


            /* New respiration flux */
            f->auto_resp =  f->gpp - f->npp;
            recalc_wb = TRUE;

            /* Now reduce LAI for down-regulated growth. */
            if (c->deciduous_model) {
                if (float_eq(s->shoot, 0.0)) {
                    s->lai = 0.0;
                } else if (s->leaf_out_days[doy] > 0.0) {
                    s->lai -= lai_inc;
                    s->lai += (f->cpleaf *
                               (p->sla * M2_AS_HA / \
                               (KG_AS_TONNES * p->cfracts)) -
                               (f->deadleaves + f->ceaten) * s->lai / s->shoot);
                } else {
                    s->lai = 0.0;
                }
            } else {
                /* update leaf area [m2 m-2] */
                if (float_eq(s->shoot, 0.0)) {
                    s->lai = 0.0;
                } else {
                    s->lai -= lai_inc;
                    s->lai += (f->cpleaf *
                               (p->sla * M2_AS_HA / \
                               (KG_AS_TONNES * p->cfracts)) -
                               (f->deadleaves + f->ceaten) * s->lai / s->shoot);
                }
            }


        }

        ntot -= f->npbranch + f->npstemimm + f->npstemmob + f->npcroot;
        ntot = MAX(0.0, ntot);

        /* allocate remaining N to flexible-ratio pools */
        f->npleaf = ntot * f->alleaf / (f->alleaf + f->alroot * p->ncrfac);
        f->nproot = ntot - f->npleaf;
    }
    return (recalc_wb);
}

//n limitation on growth
double calculate_growth_stress_limitation(params *p, state *s) {
    /* Calculate level of stress due to nitrogen or water availability */
    double nlim, current_limitation;
    double nc_opt = 0.04;

    /* N limitation based on leaf NC ratio */
    if (s->shootnc < p->nf_min) {
        nlim = 0.0;
    } else if (s->shootnc < nc_opt && s->shootnc > p->nf_min) {
        nlim = 1.0 - ((nc_opt - s->shootnc) / (nc_opt - p->nf_min));
    } else {
        nlim = 1.0;
    }

    /*
     * Limitation by nitrogen and water. Water constraint is implicit,
     * in that, water stress results in an increase of root mass,
     * which are assumed to spread horizontally within the rooting zone.
     * So in effect, building additional root mass doesnt alleviate the
     * water limitation within the model. However, it does more
     * accurately reflect an increase in root C production at a water
     * limited site. This implementation is also consistent with other
     * approaches, e.g. LPJ. In fact I dont see much evidence for models
     * that have a flexible bucket depth. Minimum constraint is limited to
     * 0.1, following Zaehle et al. 2010 (supp), eqn 18.
     */
    current_limitation = MAX(0.1, MIN(nlim, s->wtfac_root));
    return (current_limitation);
}

double alloc_goal_seek(double simulated, double target, double alloc_max,
                       double sensitivity) {

    /* Sensitivity parameter characterises how allocation fraction respond
       when the leaf:sapwood area ratio departs from the target value
       If sensitivity close to 0 then the simulated leaf:sapwood area ratio
       will closely track the target value */
    double frac = 0.5 + 0.5 * (1.0 - simulated / target) / sensitivity;

    return MAX(0.0, alloc_max * MIN(1.0, frac));
}



void precision_control(fluxes *f, state *s) {
    /* Detect very low values in state variables and force to zero to
    avoid rounding and overflow errors */

    double tolerance = 1E-10;

    /* C & N state variables */
	if (s->nsc < tolerance) {
		s->nsc = 0.0;
	}
    if (s->shoot < tolerance) {
        f->deadleaves += s->shoot;
        f->deadleafn += s->shootn;
        s->shoot = 0.0;
        s->shootn = 0.0;
    }

    if (s->branch < tolerance) {
        f->deadbranch += s->branch;
        f->deadbranchn += s->branchn;
        s->branch = 0.0;
        s->branchn = 0.0;
    }

    if (s->root < tolerance) {
        f->deadrootn += s->rootn;
        f->deadroots += s->root;
        s->root = 0.0;
        s->rootn = 0.0;
    }

    if (s->croot < tolerance) {
        f->deadcrootn += s->crootn;
        f->deadcroots += s->croot;
        s->croot = 0.0;
        s->crootn = 0.0;
    }

    /* Not setting these to zero as this just leads to errors with desert
       regrowth...instead seeding them to a small value with a CN~25. */

    if (s->stem < tolerance) {
        f->deadstems += s->stem;
        f->deadstemn += s->stemn;
        s->stem = 0.00;
        s->stemn = 0.0000;
        s->stemnimm = 0.0000;
        s->stemnmob = 0.0;
    }

    /* need separate one as this will become very small if there is no
       mobile stem N */
    if (s->stemnmob < tolerance) {
        f->deadstemn += s->stemnmob;
        s->stemnmob = 0.0;
    }

    if (s->stemnimm < tolerance) {
        f->deadstemn += s->stemnimm;
        s->stemnimm = 0.0000;
    }

    return;
}


void calculate_cn_store(control *c, fluxes *f, state *s) {
    /*
    Calculate labile C & N stores from which growth is allocated in the
    following year.
    */


    s->cstore += f->npp;
    s->nstore += f->nuptake + f->retrans;
    s->anpp += f->npp;

    /*
    double nstore_max, excess, k;
    double leaf_nc_max = 0.04;
    double CN_max = 100.0;


    if (c->alloc_model == GRASSES) {
        k = 0.3;
        nstore_max = MAX(1E-04, k * s->root * leaf_nc_max);
    } else {
        k = 0.15;
        nstore_max = MAX(1E-04, k * s->sapwood * leaf_nc_max);
    }

    s->nstore += f->nuptake + f->retrans;
    if (s->nstore > nstore_max) {
        s->nstore = nstore_max;
        f->nuptake = 0.0;
        f->retrans = 0.0;
    }

    s->cstore += f->npp;
    */
    /*
    if (s->cstore/nstore_max > CN_max) {
        excess = s->cstore - (nstore_max * CN_max);
        s->cstore = nstore_max * CN_max;
        f->auto_resp += excess;
    }*/


    return;
}


void calculate_average_alloc_fractions(fluxes *f, state *s, int days) {
    double excess;

    s->avg_alleaf /= (float) days;
    s->avg_alroot /= (float) days;
    s->avg_alcroot /= (float) days;
    s->avg_albranch /= (float) days;
    s->avg_alstem /= (float) days;

    f->alleaf = s->avg_alleaf;
    f->alroot = s->avg_alroot;
    f->alcroot = s->avg_alcroot;
    f->albranch = s->avg_albranch;
    f->alstem = s->avg_alstem;

    /*
        Because we are taking the average growing season fracs the total may
        end up being just under 1, due to rounding. So put the missing
        bit into the leaves - arbitary decision there
    */
    excess = 1.0 - f->alleaf - f->alroot - f->alcroot - f->albranch - f->alstem;
    f->alleaf += excess;

    return;
}

void allocate_stored_c_and_n(fluxes *f, params *p, state *s) {
    /*
    Allocate stored C&N. This is either down as the model is initialised
    for the first time or at the end of each year.
    */
    double ntot;

    /* ========================
       Carbon - fixed fractions
       ======================== */
    s->c_to_alloc_shoot = f->alleaf * s->cstore;
    s->c_to_alloc_root = f->alroot * s->cstore;
    s->c_to_alloc_croot = f->alcroot * s->cstore;
    s->c_to_alloc_branch = f->albranch * s->cstore;
    s->c_to_alloc_stem = f->alstem * s->cstore;

    /* =========================================================
        Nitrogen - Fixed ratios N allocation to woody components.
       ========================================================= */

    /* N flux into new ring (immobile component -> structrual components) */
    s->n_to_alloc_stemimm = s->cstore * f->alstem * p->ncwimm;

    /* N flux into new ring (mobile component -> can be retrans for new
       woody tissue) */
    s->n_to_alloc_stemmob = s->cstore * f->alstem * (p->ncwnew - p->ncwimm);
    s->n_to_alloc_branch = s->cstore * f->albranch * p->ncbnew;
    s->n_to_alloc_croot = s->cstore * f->alcroot * p->nccnew;

    /* Calculate remaining N left to allocate to leaves and roots */
    ntot = MAX(0.0, (s->nstore - s->n_to_alloc_stemimm - s->n_to_alloc_stemmob -
                     s->n_to_alloc_branch));

    /* allocate remaining N to flexible-ratio pools */
    s->n_to_alloc_shoot = (ntot * f->alleaf /
                            (f->alleaf + f->alroot * p->ncrfac));
    s->n_to_alloc_root = ntot - s->n_to_alloc_shoot;

    /*
    leaf_NC = s->n_to_alloc_shoot / s->c_to_alloc_shoot
    if leaf_NC > 0.04:
        s->n_to_alloc_shoot = s->c_to_alloc_shoot * 0.14

    s->n_to_alloc_root = ntot - s->n_to_alloc_shoot


    root_NC = s->n_to_alloc_root / s->c_to_alloc_root
    ncmaxr = 0.04 * p->ncrfac
    if root_NC > ncmaxr:
        extrar = (s->n_to_alloc_root -
                  (s->c_to_alloc_root * ncmaxr))

        s->inorgn += extrar
        s->n_to_alloc_root -= extrar
    */
    return;
}


double nitrogen_retrans(control *c, fluxes *f, params *p, state *s,
                        double fdecay, double rdecay, int doy) {
    /* Nitrogen retranslocated from senesced plant matter.
    Constant rate of n translocated from mobile pool

    Parameters:
    -----------
    fdecay : float
        foliage decay rate
    rdecay : float
        fine root decay rate

    Returns:
    --------
    N retrans : float
        N retranslocated plant matter

    */
    double leafretransn, rootretransn, crootretransn, branchretransn,
           stemretransn;

    if (c->deciduous_model) {
        leafretransn = p->fretrans * f->lnrate * s->remaining_days[doy];
    } else {
        leafretransn = p->fretrans * fdecay * s->shootn;
    }

    rootretransn = p->rretrans * rdecay * s->rootn;
    crootretransn = p->cretrans * p->crdecay * s->crootn;
    branchretransn = p->bretrans * p->bdecay * s->branchn;
    stemretransn = (p->wretrans * p->wdecay * s->stemnmob + p->retransmob *
                    s->stemnmob);

    /* store for NCEAS output */
    f->leafretransn = leafretransn;

    return (leafretransn + rootretransn + crootretransn + branchretransn +
            stemretransn);
}


double calculate_nuptake(control *c, params *p, state *s) {
    /* N uptake depends on the rate at which soil mineral N is made
    available to the plants.

    Returns:
    --------
    nuptake : float
        N uptake

    References:
    -----------
    * Dewar and McMurtrie, 1996, Tree Physiology, 16, 161-171.
    * Raich et al. 1991, Ecological Applications, 1, 399-429.

    */
    double nuptake, U0, Kr;

    if (c->nuptake_model == 0) {
        /* Constant N uptake */
        nuptake = p->nuptakez;
    } else if (c->nuptake_model == 1) {
        /* evaluate nuptake : proportional to dynamic inorganic N pool */
        nuptake = p->rateuptake * s->inorgn;
    } else if (c->nuptake_model == 2) {
        /* N uptake is a saturating function on root biomass following
           Dewar and McMurtrie, 1996. */

        /* supply rate of available mineral N */
        U0 = p->rateuptake * s->inorgn;
        Kr = p->kr;
        nuptake = MAX(U0 * s->root / (s->root + Kr), 0.0);

        /* Make minimum uptake rate supply rate for deciduous_model cases
           otherwise it is possible when growing from scratch we don't have
           enough root mass to obtain N at the annual time step
           I don't see an obvious better solution?
        if c->deciduous_model:
            nuptake = max(U0 * s->root / (s->root + Kr), U0) */
    } else {
        fprintf(stderr, "Unknown N uptake option\n");
        exit(EXIT_FAILURE);
    }

    return (nuptake);
}


void initialise_roots(fluxes *f, params *p, state *s) {
    /* Set up all the rooting arrays for use with the hydraulics assumptions */
    int    i;
    double thick;

    // Using CABLE depths, but spread over 2 m.
    //double cable_thickness[7] = {0.01, 0.025, 0.067, 0.178, 0.472, 1.248, 2.0};

    s->thickness = malloc(p->core * sizeof(double));
    if (s->thickness == NULL) {
        fprintf(stderr, "malloc failed allocating thickness\n");
        exit(EXIT_FAILURE);
    }

    /* root mass is g biomass, i.e. ~twice the C content */
    s->root_mass = malloc(p->core * sizeof(double));
    if (s->root_mass == NULL) {
        fprintf(stderr, "malloc failed allocating root_mass\n");
        exit(EXIT_FAILURE);
    }

    s->root_length = malloc(p->core * sizeof(double));
    if (s->root_length == NULL) {
        fprintf(stderr, "malloc failed allocating root_length\n");
        exit(EXIT_FAILURE);
    }

    s->layer_depth = malloc(p->core * sizeof(double));
    if (s->layer_depth == NULL) {
        fprintf(stderr, "malloc failed allocating layer_depth\n");
        exit(EXIT_FAILURE);
    }

    // force a thin top layer = 0.1
    thick = 0.1;
    s->layer_depth[0] = thick;
    s->thickness[0] = thick;
    for (i = 1; i < p->core; i++) {
        thick += p->layer_thickness;
        s->layer_depth[i] = thick;
        s->thickness[i] = p->layer_thickness;

        /* made up initalisation, following SPA, get replaced second timestep */
        s->root_mass[i] = 0.1;
        s->root_length[i] = 0.1;
        s->rooted_layers = p->core;
    }

    return;
}

void update_roots(control *c, params *p, state *s) {
    /*
        Given the amount of roots grown by GDAY predict the assoicated rooting
        distribution accross soil layers
        - These assumptions come from Mat's SPA model.

        TODO: implement CABLE version.
    */
    int    i;
    double C_2_BIOMASS = 2.0;
    double min_biomass;
    double root_biomass, root_cross_sec_area, root_depth, root_reach, mult;
    double surf_biomass, prev, curr, slope, cumulative_depth;
    double x1 = 0.1;        /* lower bound for brent search */
    double x2 = 10.0;       /* upper bound for brent search */
    double tol = 0.0001;    /* tolerance for brent search */
    double fine_root, fine_root_min;

    min_biomass = 20.0;  // g root biomss
    fine_root = s->root * TONNES_HA_2_G_M2 * C_2_BIOMASS;
    root_biomass = MAX(min_biomass, fine_root);
    //root_biomass = MAX(min_biomass,  305.0 * C_2_BIOMASS);

    root_cross_sec_area = M_PI * p->root_radius * p->root_radius;   /* (m2) */
    root_depth = p->max_depth * root_biomass / (p->root_k + root_biomass);

    s->rooted_layers = 0;
    for (i = 0; i < p->core; i++) {
        if (s->layer_depth[i] > root_depth) {
            s->rooted_layers = i;
            break;
        }
    }

    /* how for into the soil do the reach extend? */
    root_reach = s->layer_depth[s->rooted_layers];

    /* Enforce 50 % of root mass in the top 1/4 of the rooted layers. */
    mult = MIN(1.0 / s->thickness[0], \
               MAX(2.0, 11.0 * exp(-0.006 * root_biomass)));

    /*
    ** assume surface root density is related to total root biomass by a
    ** scalar dependent on root biomass
    */
    surf_biomass = root_biomass * mult;

    if (s->rooted_layers > 1) {
        /*
        ** determine slope of root distribution given rooting depth
        ** and ratio of root mass to surface root density
        */
        slope = zbrent(&calc_root_dist, x1, x2, tol, root_biomass,
                       surf_biomass, s->rooted_layers, s->thickness[0],
                       root_reach);

        prev = 1.0 / slope;
        cumulative_depth = 0.0;
        for (i = 0; i <= s->rooted_layers; i++) {
            cumulative_depth += s->thickness[i];
            curr = 1.0 / slope * exp(-slope * cumulative_depth);
            s->root_mass[i] = (prev - curr) * surf_biomass;

            /* (m m-3 soil) */
            s->root_length[i] = s->root_mass[i] / (p->root_density * \
                                                   root_cross_sec_area);
            prev = curr;
        }
    } else {
        s->root_mass[0] = root_biomass;
    }

    return;
}

double calc_root_dist(double slope, double root_biomass, double surf_biomass,
                      double rooted_layers, double top_lyr_thickness,
                      double root_reach) {
    /*
        This function is used in the in the zbrent numerical algorithm to
        figure out the slope of the rooting distribution for a given depth
    */
    double one, two, arg1, arg2;

    one = (1.0 - exp(-slope * rooted_layers * top_lyr_thickness)) / slope;
    two = root_biomass / surf_biomass;
    arg1 = (1.0 - exp(-slope * root_reach)) / slope;
    arg2 = root_biomass / surf_biomass;

    return (arg1 - arg2);
}
