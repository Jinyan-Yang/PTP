/* ============================================================================
* Calculate C and N litter production
*
* Litter production for each pool is assumed to be proportional to biomass
* pool size.
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
#include "litter_production.h"

void calculate_litterfall(control *c, fluxes *f, fast_spinup *fs,
                          params *p, state *s, int doy, double *fdecay,
                          double *rdecay) {

    double  ncflit, ncrlit;

    /* Leaf/root litter rates are higher during dry periods and therefore is
    dependent on soil water content */
    *fdecay = decay_in_dry_soils(p->fdecay, p->fdecaydry, p, s);
    *rdecay = decay_in_dry_soils(p->rdecay, p->rdecaydry, p, s);

    /* litter N:C ratios, roots and shoot */
    ncflit = s->shootnc * (1.0 - p->fretrans);
    ncrlit = s->rootnc * (1.0 - p->rretrans);

    /* C litter production */
    f->deadroots = *rdecay * s->root;
    f->deadcroots = p->crdecay * s->croot;
    f->deadstems = p->wdecay * s->stem;
    f->deadbranch = p->bdecay * s->branch;
    f->deadsapwood = (p->wdecay + p->sapturnover) * s->sapwood;


    //if (c->deciduous_model)
    //    f->deadleaves = f->lrate * s->remaining_days[doy];
    //else
        //f->deadleaves = *fdecay * s->shoot;
	
	//f->deadleaves = (p->fdecay) * s->shoot;
    f->deadleaves = p->fdecay* pow((1 - s->wtfac_topsoil), p->q_s)* s->shoot;

    if (c->spinup_method == SAS) {
        if (c->deciduous_model)
            fs->loss[LF] += f->lrate;
        else
            fs->loss[LF] += *fdecay;

        fs->loss[LR] += *rdecay;
        fs->loss[LCR] += p->crdecay;
        fs->loss[LB] += p->bdecay;
        fs->loss[LW] += p->wdecay;
    }

    /* N litter production */
    f->deadleafn = f->deadleaves * ncflit;

    /* Assuming fraction is retranslocated before senescence, i.e. a fracion
       of nutrients is stored within the plant */
    f->deadrootn = f->deadroots * ncrlit;
    f->deadcrootn = p->crdecay * s->crootn * (1.0 - p->cretrans);
    f->deadbranchn = p->bdecay * s->branchn * (1.0 - p->bretrans);

    /* N in stemwood litter - only mobile n is retranslocated */
    f->deadstemn = p->wdecay * (s->stemnimm + s->stemnmob * \
                    (1.0 - p->wretrans));

    /* Animal grazing? */

    /* Daily... */
//    if (c->grazing == 1) {
//        daily_grazing_calc(*fdecay, p, f, s);
//
//    /* annually */
//    } else if (c->grazing == 2 && p->disturbance_doy == doy) {
//        annual_grazing_calc(p, f, s);
//
//    /* no grazing */
//    } else {
//        f->ceaten = 0.0;
//        f->neaten = 0.0;
//    }
    return;

}
void calculate_harvest(fluxes *f, params *p, state *s, int doy, int year) {
    /*
    // do grazing/harvest in a cheating way
    const char* tmp_year = p->year_harvest;
    const char* tmp_doy = p->doy_harvest;
    char resultYear[4];
    char resultDOY[3];
    // assign current yr and doy date as char
    sprintf(resultYear, "%d", year);
    sprintf(resultDOY, "%d", doy);

    // check if current doy and yr is in the input char
    char* ret_yr;
    char* ret_doy;
    ret_yr = strstr(tmp_year, resultYear);
    ret_doy = strstr(tmp_doy, resultDOY);
    // reduce aboveground biomass when there is an harvest
    // currently the fraction of reduction is a arbitory
    if (ret_yr && ret_doy) {

       // f->ceaten = 0.8 * s->shoot;
        //f->neaten = 0.8 * s->shootn;

        f->ceaten = 0.0;
        f->neaten = 0.0;

        //set harect c
        s->shoot = 0.2 * s->shoot;
    }
    else {
        f->ceaten = 0;
        f->neaten = 0;
    }
    */

    // do grazing/harvest in a cheating way
    char tmp_year[255];
    strcpy(tmp_year, p->year_harvest);
  //  char tmp_doy[255];
   // strcpy(tmp_doy, p->doy_harvest);

    char resultYear[9];
    char resultDOY[4];
    // assign current yr and doy date as char
    sprintf(resultYear, "%d", year);
    sprintf(resultDOY, "%d", doy);

    strcat(resultYear, resultDOY);
    strcat(resultYear, ",");
    // check if current doy and yr is in the input char
    char *ret_yr;

    ret_yr = strstr(tmp_year, resultYear);
   // ret_doy = strstr(tmp_doy, resultDOY);
    // reduce aboveground biomass when there is an harvest
    // currently the fraction of reduction is a arbitory
    if (ret_yr) {

        // f->ceaten = 0.8 * s->shoot;
         //f->neaten = 0.8 * s->shootn;

        f->ceaten = 0.0;
        f->neaten = 0.0;

        //set harect c
        s->shoot = 0.2 * s->shoot;
    }
    else {
        f->ceaten = 0;
        f->neaten = 0;
    }


    return(0);
}
void daily_grazing_calc(double fdecay, params *p, fluxes *f, state *s) {
    /* daily grass grazing...

    Parameters:
    -----------
    fdecay : float
        foliage decay rate

    Returns:
    --------
    ceaten : float
        C consumed by grazers [tonnes C/ha/day]
    neaten : float
        N consumed by grazers [tonnes C/ha/day]
    */
    f->ceaten = fdecay * p->fracteaten / (1.0 - p->fracteaten) * s->shoot;
    f->neaten = fdecay * p->fracteaten / (1.0 - p->fracteaten) * s->shootn;

    return;
}

void annual_grazing_calc(params *p, fluxes *f, state *s) {
    /* Annual grass grazing...single one off event


    Returns:
    --------
    ceaten : float
        C consumed by grazers [tonnes C/ha/day]
    neaten : float
        N consumed by grazers [tonnes C/ha/day]
    */
    f->ceaten = s->shoot * p->fracteaten;
    f->neaten = s->shootn * p->fracteaten;

    return;
}

float decay_in_dry_soils(double decay_rate, double decay_rate_dry, params *p,
                         state *s) {
    /* Decay rates (e.g. leaf litterfall) can increase in dry soil, adjust
    decay param. This is based on field measurements by F. J. Hingston
    (unpublished) cited in Corbeels.

    Parameters:
    -----------
    decay_rate : float
        default model parameter decay rate [tonnes C/ha/day]
    decay_rate_dry : float
        default model parameter dry deacy rate [tonnes C/ha/day]

    Returns:
    --------
    decay_rate : float
        adjusted deacy rate if the soil is dry [tonnes C/ha/day]

    Reference:
    ----------
    Corbeels et al. (2005) Ecological Modelling, 187, 449-474.

    */
    /* turn into fraction... */

	// changed to a simpler decay function based on sw
    double smc_root, new_decay_rate;
    //smc_root = s->pawater_root / p->wcapac_root;

    //new_decay_rate = (decay_rate_dry - (decay_rate_dry - decay_rate) *
    //                 (smc_root - p->watdecaydry) /
    //                 (p->watdecaywet - p->watdecaydry));

    //if (new_decay_rate < decay_rate)
    //    new_decay_rate = decay_rate;

    //if (new_decay_rate > decay_rate_dry)
    //    new_decay_rate = decay_rate_dry;

	new_decay_rate = decay_rate_dry * (1-s->wtfac_root+0.1);

    return new_decay_rate;
}
