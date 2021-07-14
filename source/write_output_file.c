/* ============================================================================
* Print output file (ascii/binary)
*
*
*
* NOTES:
*   Currently I have only implemented the ascii version
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   25.02.2015
*
* =========================================================================== */
#include "write_output_file.h"


void open_output_file(control *c, char *fname, FILE **fp) {
    *fp = fopen(fname, "w");
    if (*fp == NULL)
        prog_error("Error opening output file for write on line", __LINE__);
}

void write_output_subdaily_header(control *c, FILE **fp) {
    /*
        Write 30 min fluxes headers to an output CSV file. This is very basic
        for now...
    */

    ///* Git version */
    //fprintf(*fp, "#Git_revision_code:%s\n", c->git_code_ver);

    /* time stuff */
    fprintf(*fp, "year,doy,hod,");

    /*
    ** Canopy stuff...
    */
    fprintf(*fp, "an_canopy,rd_canopy,gsc_canopy,");
    fprintf(*fp, "apar_canopy,trans_canopy,tleaf\n");
    return;
}

void write_output_header(control *c, FILE **fp) {
    /*
        Write daily state and fluxes headers to an output CSV file. Note we
        are not writing anything useful like units as there is a wrapper
        script to translate the outputs to a nice CSV file with input met
        data, units and nice header information.
    */
    int ncols = 86;
    int nrows = c->num_days;

    ///* Git version */
    //fprintf(*fp, "#Git_revision_code:%s\n", c->git_code_ver);

    /* time stuff */
    fprintf(*fp, "year,doy,");

    /*
    ** STATE
    */

    /* water*/
    fprintf(*fp, "wtfac_root,wtfac_topsoil,pawater_root,");

    /* plant */
    fprintf(*fp, "nsc,shoot,lai,branch,stem,root,croot,shootn,branchn,stemn,");
    fprintf(*fp, "rootn,crootn,cstore,nstore,");

    /* belowground */
    fprintf(*fp, "soilc,soiln,inorgn,litterc,littercag,littercbg,");
    fprintf(*fp, "litternag,litternbg,activesoil,slowsoil,");
    fprintf(*fp, "passivesoil,activesoiln,slowsoiln,passivesoiln,");

    /*
    ** FLUXES
    */

    /* water */
    fprintf(*fp, "et,transpiration,soil_evap,canopy_evap,runoff,");
    fprintf(*fp, "gs_mol_m2_sec,ga_mol_m2_sec,");

    /* litter */
    fprintf(*fp, "deadleaves,deadbranch,deadstems,deadroots,deadcroots,");
    fprintf(*fp, "deadleafn,deadbranchn,deadstemn,deadrootn,deadcrootn,");

    /* C fluxes */
    fprintf(*fp, "nep,gpp,npp,hetero_resp,auto_resp,apar,");

    /* C & N growth */
    fprintf(*fp, "cpleaf,cpbranch,cpstem,cproot,cpcroot,");
    fprintf(*fp, "npleaf,npbranch,npstemimm,npstemmob,nproot,npcroot,");

    /* N stuff */
    fprintf(*fp, "nuptake,ngross,nmineralisation,nloss,");

    /* traceability stuff */
    fprintf(*fp, "tfac_soil_decomp,c_into_active,c_into_slow,");
    fprintf(*fp, "c_into_passive,active_to_slow,active_to_passive,");
    fprintf(*fp, "slow_to_active,slow_to_passive,passive_to_active,");
    fprintf(*fp, "co2_rel_from_surf_struct_litter,");
    fprintf(*fp, "co2_rel_from_soil_struct_litter,");
    fprintf(*fp, "co2_rel_from_surf_metab_litter,");
    fprintf(*fp, "co2_rel_from_soil_metab_litter,");
    fprintf(*fp, "co2_rel_from_active_pool,");
    fprintf(*fp, "co2_rel_from_slow_pool,");
    fprintf(*fp, "co2_rel_from_passive_pool,");

    /* extra priming stuff */
    fprintf(*fp, "root_exc,");
    fprintf(*fp, "root_exn,");
    fprintf(*fp, "co2_released_exud,");
    fprintf(*fp, "factive,");
    fprintf(*fp, "rtslow,");
    fprintf(*fp, "rexc_cue,");



    /* Misc */
    fprintf(*fp, "predawn_swp,");
    fprintf(*fp, "midday_lwp,");
    fprintf(*fp, "midday_xwp,");
    fprintf(*fp, "leafretransn,");
    fprintf(*fp, "dead_year,");
    fprintf(*fp, "dead_doy,");
    fprintf(*fp, "theta0,");
    fprintf(*fp, "theta1,");
    fprintf(*fp, "theta2,");
    fprintf(*fp, "theta3,");
    fprintf(*fp, "theta4,");
    fprintf(*fp, "theta5,");
    fprintf(*fp, "theta6,");
    fprintf(*fp, "theta7,");
    fprintf(*fp, "theta8,");
    fprintf(*fp, "theta9,");
    fprintf(*fp, "theta10,");
    fprintf(*fp, "theta11,");
    fprintf(*fp, "theta12,");
    fprintf(*fp, "theta13,");
    fprintf(*fp, "theta14,");
    fprintf(*fp, "theta15,");
    fprintf(*fp, "theta16,");
    fprintf(*fp, "theta17,");
    fprintf(*fp, "theta18,");
    fprintf(*fp, "theta19,");
    fprintf(*fp, "theta20\n");

    if (c->output_ascii == FALSE) {
        fprintf(*fp, "nrows=%d\n", nrows);
        fprintf(*fp, "ncols=%d\n", ncols);
    }
    return;
}

void write_subdaily_outputs_ascii(control *c, canopy_wk *cw, double year,
                                  double doy, int hod) {
    /*
        Write sub-daily canopy fluxes - very basic for now
    */

    /* time stuff */
    fprintf(c->ofp_sd, "%.10f,%.10f,%.10f,", year, doy, (double)hod);

    /* Canopy stuff */
    fprintf(c->ofp_sd, "%.10f,%.10f,%.10f,",
                       cw->an_canopy, cw->rd_canopy, cw->gsc_canopy);
    fprintf(c->ofp_sd, "%.10f,%.10f,%.10f\n",
                       cw->apar_canopy, cw->trans_canopy, cw->tleaf_new);

    return;
}

void write_daily_outputs_ascii(control *c, canopy_wk *cw, fluxes *f, state *s,
                               int year, int doy) {
    /*
        Write daily state and fluxes headers to an output CSV file. Note we
        are not writing anything useful like units as there is a wrapper
        script to translate the outputs to a nice CSV file with input met
        data, units and nice header information.
    */


    /* time stuff */
    fprintf(c->ofp, "%.10f,%.10f,", (double)year, (double)doy);

    /*
    ** STATE

    */

    /* water*/
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
            s->wtfac_root,s->wtfac_topsoil,s->pawater_root);

    /* plant */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
		            s->nsc, s->shoot,s->lai,s->branch,s->stem,s->root,s->croot,
                    s->shootn,s->branchn,s->stemn);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    s->rootn,s->crootn,s->cstore,s->nstore);

    /* belowground */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
                    s->soilc,s->soiln,s->inorgn,s->litterc,s->littercag,
                    s->littercbg);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
                    s->litternag,s->litternbg,s->activesoil,s->slowsoil,
                    s->passivesoil,s->activesoiln,s->slowsoiln,
                    s->passivesoiln);
    /*
    ** FLUXES
    */

    /* water */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    f->et,f->transpiration,f->soil_evap,f->canopy_evap);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
                    f->runoff,f->gs_mol_m2_sec,f->ga_mol_m2_sec);

    /* litter */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->deadleaves,f->deadbranch,f->deadstems,f->deadroots,
                    f->deadcroots);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->deadleafn,f->deadbranchn,f->deadstemn,f->deadrootn,
                    f->deadcrootn);

    /* C fluxes */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->nep,f->gpp,f->npp,f->hetero_resp,f->auto_resp,
                    f->apar);

    /* C & N growth */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->cpleaf,f->cpbranch,f->cpstem,f->cproot,f->cpcroot);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,",
                    f->npleaf,f->npbranch,f->npstemimm,f->npstemmob,f->nproot,
                    f->npcroot);

    /* N stuff */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,%.10f,",
                    f->nuptake,f->ngross,f->nmineralisation,f->nloss);


    /* traceability stuff */
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
                    f->tfac_soil_decomp,f->c_into_active,f->c_into_slow);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
                    f->c_into_passive,f->active_to_slow,f->active_to_passive);
    fprintf(c->ofp, "%.10f,%.10f,%.10f,",
                    f->slow_to_active,f->slow_to_passive,f->passive_to_active);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_surf_struct_litter);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_soil_struct_litter);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_surf_metab_litter);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_soil_metab_litter);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_active_pool);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_slow_pool);
    fprintf(c->ofp, "%.10f,", f->co2_rel_from_passive_pool);

    /* extra priming stuff */
    fprintf(c->ofp, "%.10f,", f->root_exc);
    fprintf(c->ofp, "%.10f,", f->root_exn);
    fprintf(c->ofp, "%.10f,", f->co2_released_exud);
    fprintf(c->ofp, "%.10f,", f->factive);
    fprintf(c->ofp, "%.10f,", f->rtslow);
    fprintf(c->ofp, "%.10f,", f->rexc_cue);

    /* Misc */
    fprintf(c->ofp, "%.10f,", s->predawn_swp);
    fprintf(c->ofp, "%.10f,", s->midday_lwp);
    fprintf(c->ofp, "%.10f,", s->midday_xwp);
    fprintf(c->ofp, "%.10f,", f->leafretransn);
    fprintf(c->ofp, "%.10f,", cw->death_year);
    fprintf(c->ofp, "%.10f,", cw->death_doy);

    if (c->water_balance == HYDRAULICS) {
        fprintf(c->ofp, "%.10f,", s->water_frac[0]);
        fprintf(c->ofp, "%.10f,", s->water_frac[1]);
        fprintf(c->ofp, "%.10f,", s->water_frac[2]);
        fprintf(c->ofp, "%.10f,", s->water_frac[3]);
        fprintf(c->ofp, "%.10f,", s->water_frac[4]);
        fprintf(c->ofp, "%.10f,", s->water_frac[5]);
        fprintf(c->ofp, "%.10f,", s->water_frac[6]);
        fprintf(c->ofp, "%.10f,", s->water_frac[7]);
        fprintf(c->ofp, "%.10f,", s->water_frac[8]);
        fprintf(c->ofp, "%.10f,", s->water_frac[9]);
        fprintf(c->ofp, "%.10f,", s->water_frac[10]);
        fprintf(c->ofp, "%.10f,", s->water_frac[11]);
        fprintf(c->ofp, "%.10f,", s->water_frac[12]);
        fprintf(c->ofp, "%.10f,", s->water_frac[13]);
        fprintf(c->ofp, "%.10f,", s->water_frac[14]);
        fprintf(c->ofp, "%.10f,", s->water_frac[15]);
        fprintf(c->ofp, "%.10f,", s->water_frac[16]);
        fprintf(c->ofp, "%.10f,", s->water_frac[17]);
        fprintf(c->ofp, "%.10f,", s->water_frac[18]);
        fprintf(c->ofp, "%.10f,", s->water_frac[19]);
        fprintf(c->ofp, "%.10f\n", s->water_frac[20]);
    } else {
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f,", -999.9);
        fprintf(c->ofp, "%.10f\n", -999.9);
    }

    return;
}

void write_daily_outputs_binary(control *c, fluxes *f, state *s, int year,
                                int doy) {
    /*
        Write daily state and fluxes headers to an output CSV file. Note we
        are not writing anything useful like units as there is a wrapper
        script to translate the outputs to a nice CSV file with input met
        data, units and nice header information.
    */
    double temp;

    /* time stuff */
    temp = (double)year;
    fwrite(&temp, sizeof(double), 1, c->ofp);
    temp = (double)doy;
    fwrite(&temp, sizeof(double), 1, c->ofp);


    /* plant */
    fwrite(&(s->shoot), sizeof(double), 1, c->ofp);
    fwrite(&(s->lai), sizeof(double), 1, c->ofp);
    fwrite(&(s->branch), sizeof(double), 1, c->ofp);
    fwrite(&(s->stem), sizeof(double), 1, c->ofp);
    fwrite(&(s->root), sizeof(double), 1, c->ofp);

    /* water */
    fwrite(&(s->wtfac_root), sizeof(double), 1, c->ofp);
    fwrite(&(s->pawater_root), sizeof(double), 1, c->ofp);
    fwrite(&(f->transpiration), sizeof(double), 1, c->ofp);
    fwrite(&(f->soil_evap), sizeof(double), 1, c->ofp);
    fwrite(&(f->canopy_evap), sizeof(double), 1, c->ofp);
    fwrite(&(f->runoff), sizeof(double), 1, c->ofp);

    /* C fluxes */
    fwrite(&(f->npp), sizeof(double), 1, c->ofp);

    return;
}


int write_final_state(control *c, params *p, state *s)
{
    /*
    Write the final state to the input param file so we can easily restart
    the model. This function copies the input param file with the exception
    of anything in the git hash and the state which it replaces with the updated
    stuff.

    */

    char line[STRING_LENGTH];
    char saved_line[STRING_LENGTH];
    char section[STRING_LENGTH] = "";
    char prev_name[STRING_LENGTH] = "";
    char *start;
    char *end;
    char *name;
    char *value;

    int error = 0;
    int line_number = 0;
    int match = FALSE;

    while (fgets(line, sizeof(line), c->ifp) != NULL) {
        strcpy(saved_line, line);
        line_number++;
        start = lskip(rstrip(line));
        if (*start == ';' || *start == '#') {
            /* Per Python ConfigParser, allow '#' comments at start of line */
        }
        else if (*start == '[') {
            /* A "[section]" line */
            end = find_char_or_comment(start + 1, ']');
            if (*end == ']') {
                *end = '\0';
                strncpy0(section, start + 1, sizeof(section));
                *prev_name = '\0';

            }
            else if (!error) {
                /* No ']' found on section line */
                error = line_number;

            }
        }
        else if (*start && *start != ';') {
            /* Not a comment, must be a name[=:]value pair */
            end = find_char_or_comment(start, '=');
            if (*end != '=') {
                end = find_char_or_comment(start, ':');
            }
            if (*end == '=' || *end == ':') {
                *end = '\0';
                name = rstrip(start);
                value = lskip(end + 1);
                end = find_char_or_comment(value, '\0');
                if (*end == ';')
                    *end = '\0';
                rstrip(value);

                /* Valid name[=:]value pair found, call handler */
                strncpy0(prev_name, name, sizeof(prev_name));

                if (!ohandler(section, name, value, c, p, s, &match) && !error)
                    error = line_number;
            }
            else if (!error) {
                /* No '=' or ':' found on name[=:]value line */
                error = line_number;
                break;
            }
        }
        if (match == FALSE)
            fprintf(c->ofp, "%s", saved_line);
        else
            match = FALSE; /* reset match flag */
    }
    return error;

}


int ohandler(char *section, char *name, char *value, control *c, params *p,
             state *s, int *match)
{
    /*
    Search for matches of the git and state values and where found write the
    current state values to the output parameter file.

    - also added previous ncd as this potential can be changed internally
    */

    #define MATCH(s, n) strcasecmp(section, s) == 0 && strcasecmp(name, n) == 0

    /*
    ** GIT
    */
    //if (MATCH("git", "git_hash")) {
    //    fprintf(c->ofp, "git_hash = %s\n", c->git_code_ver);
    //    *match = TRUE;
    //}

    /*
    ** PARAMS
    */
    if (MATCH("params", "previous_ncd")) {
        fprintf(c->ofp, "previous_ncd = %.10f\n", p->previous_ncd);
        *match = TRUE;
    }

    /*
    ** STATE
    */

    if (MATCH("state", "activesoil")) {
        fprintf(c->ofp, "activesoil = %.10f\n", s->activesoil);
        *match = TRUE;
    } else if (MATCH("state", "activesoiln")) {
        fprintf(c->ofp, "activesoiln = %.10f\n", s->activesoiln);
        *match = TRUE;
    } else if (MATCH("state", "age")) {
        fprintf(c->ofp, "age = %.10f\n", s->age);
        *match = TRUE;
    } else if (MATCH("state", "avg_albranch")) {
        fprintf(c->ofp, "avg_albranch = %.10f\n", s->avg_albranch);
        *match = TRUE;
    } else if (MATCH("state", "avg_alcroot")) {
        fprintf(c->ofp, "avg_alcroot = %.10f\n", s->avg_alcroot);
        *match = TRUE;
    } else if (MATCH("state", "avg_alleaf")) {
        fprintf(c->ofp, "avg_alleaf = %.10f\n", s->avg_alleaf);
        *match = TRUE;
    } else if (MATCH("state", "avg_alroot")) {
        fprintf(c->ofp, "avg_alroot = %.10f\n", s->avg_alroot);
        *match = TRUE;
    } else if (MATCH("state", "avg_alstem")) {
        fprintf(c->ofp, "avg_alstem = %.10f\n", s->avg_alstem);
        *match = TRUE;
    } else if (MATCH("state", "branch")) {
        fprintf(c->ofp, "branch = %.10f\n", s->branch);
        *match = TRUE;
    } else if (MATCH("state", "branchn")) {
        fprintf(c->ofp, "branchn = %.10f\n", s->branchn);
        *match = TRUE;
    } else if (MATCH("state", "canht")) {
        fprintf(c->ofp, "canht = %.10f\n", s->canht);
        *match = TRUE;
    } else if (MATCH("state", "croot")) {
        fprintf(c->ofp, "croot = %.10f\n", s->croot);
        *match = TRUE;
    } else if (MATCH("state", "crootn")) {
        fprintf(c->ofp, "crootn = %.10f\n", s->crootn);
        *match = TRUE;
    } else if (MATCH("state", "cstore")) {
        fprintf(c->ofp, "cstore = %.10f\n", s->cstore);
        *match = TRUE;
    } else if (MATCH("state", "inorgn")) {
        fprintf(c->ofp, "inorgn = %.10f\n", s->inorgn);
        *match = TRUE;
    } else if (MATCH("state", "lai")) {
        fprintf(c->ofp, "lai = %.10f\n", s->lai);
        *match = TRUE;
    } else if (MATCH("state", "metabsoil")) {
        fprintf(c->ofp, "metabsoil = %.10f\n", s->metabsoil);
        *match = TRUE;
    } else if (MATCH("state", "metabsoiln")) {
        fprintf(c->ofp, "metabsoiln = %.10f\n", s->metabsoiln);
        *match = TRUE;
    } else if (MATCH("state", "metabsurf")) {
        fprintf(c->ofp, "metabsurf = %.10f\n", s->metabsurf);
        *match = TRUE;
    } else if (MATCH("state", "metabsurfn")) {
        fprintf(c->ofp, "metabsurfn = %.10f\n", s->metabsurfn);
        *match = TRUE;
    } else if (MATCH("state", "nstore")) {
        fprintf(c->ofp, "nstore = %.10f\n", s->nstore);
        *match = TRUE;
    } else if (MATCH("state", "passivesoil")) {
        fprintf(c->ofp, "passivesoil = %.10f\n", s->passivesoil);
        *match = TRUE;
    } else if (MATCH("state", "passivesoiln")) {
        fprintf(c->ofp, "passivesoiln = %.10f\n", s->passivesoiln);
        *match = TRUE;
    } else if (MATCH("state", "pawater_root")) {
        fprintf(c->ofp, "pawater_root = %.10f\n", s->pawater_root);
        *match = TRUE;
    } else if (MATCH("state", "pawater_topsoil")) {
        fprintf(c->ofp, "pawater_topsoil = %.10f\n", s->pawater_topsoil);
        *match = TRUE;
    } else if (MATCH("state", "prev_sma")) {
        fprintf(c->ofp, "prev_sma = %.10f\n", s->prev_sma);
        *match = TRUE;
    } else if (MATCH("state", "root")) {
        fprintf(c->ofp, "root = %.10f\n", s->root);
        *match = TRUE;
    } else if (MATCH("state", "root_depth")) {
        fprintf(c->ofp, "root_depth = %.10f\n", s->root_depth);
        *match = TRUE;
    } else if (MATCH("state", "rootn")) {
        fprintf(c->ofp, "rootn = %.10f\n", s->rootn);
        *match = TRUE;
    } else if (MATCH("state", "sapwood")) {
        fprintf(c->ofp, "sapwood = %.10f\n", s->sapwood);
        *match = TRUE;
    } else if (MATCH("state", "shoot")) {
        fprintf(c->ofp, "shoot = %.10f\n", s->shoot);
        *match = TRUE;
    } else if (MATCH("state", "shootn")) {
        fprintf(c->ofp, "shootn = %.10f\n", s->shootn);
        *match = TRUE;
    } else if (MATCH("state", "sla")) {
        fprintf(c->ofp, "sla = %.10f\n", s->sla);
        *match = TRUE;
    } else if (MATCH("state", "slowsoil")) {
        fprintf(c->ofp, "slowsoil = %.10f\n", s->slowsoil);
        *match = TRUE;
    } else if (MATCH("state", "slowsoiln")) {
        fprintf(c->ofp, "slowsoiln = %.10f\n", s->slowsoiln);
        *match = TRUE;
    } else if (MATCH("state", "stem")) {
        fprintf(c->ofp, "stem = %.10f\n", s->stem);
        *match = TRUE;
    } else if (MATCH("state", "stemn")) {
        fprintf(c->ofp, "stemn = %.10f\n", s->stemn);
        *match = TRUE;
    } else if (MATCH("state", "stemnimm")) {
        fprintf(c->ofp, "stemnimm = %.10f\n", s->stemnimm);
        *match = TRUE;
    } else if (MATCH("state", "stemnmob")) {
        fprintf(c->ofp, "stemnmob = %.10f\n", s->stemnmob);
        *match = TRUE;
    } else if (MATCH("state", "structsoil")) {
        fprintf(c->ofp, "structsoil = %.10f\n", s->structsoil);
        *match = TRUE;
    } else if (MATCH("state", "structsoiln")) {
        fprintf(c->ofp, "structsoiln = %.10f\n", s->structsoiln);
        *match = TRUE;
    } else if (MATCH("state", "structsurf")) {
        fprintf(c->ofp, "structsurf = %.10f\n", s->structsurf);
        *match = TRUE;
    } else if (MATCH("state", "structsurfn")) {
        fprintf(c->ofp, "structsurfn = %.10f\n", s->structsurfn);
        *match = TRUE;
    }

    return (1);
}
