#ifndef PLANT_GROWTH_H
#define PLANT_GROWTH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


#include "gday.h"
#include "constants.h"
#include "water_balance.h"
#include "utilities.h"
#include "photosynthesis.h"
#include "optimal_root_model.h"
#include "canopy.h"

/* C stuff */
void    calc_day_growth(canopy_wk *, control *, fluxes *, fast_spinup *,
                        met_arrays *ma, met *, nrutil *, params *, state *,
                        double, int, double, double);
void    carbon_allocation(control *, fluxes *, params *, state *,
                                                     double, int);
void    calc_carbon_allocation_fracs(control *c, fluxes *, fast_spinup *,
                                     params *, state *, met_arrays* ma,double);
double  alloc_goal_seek(double, double, double, double);
void    update_plant_state(control *, fluxes *, params *, state *,
                                                        double, double, int);
void    precision_control(fluxes *, state *);
void    calculate_cn_store(control *, fluxes *, state *);
void    calculate_average_alloc_fractions(fluxes *, state *, int );
void    allocate_stored_c_and_n(fluxes *f, params *p, state *s);
void    carbon_daily_production(control *, fluxes *, met *m, params *, state *,
                                double);
void    calculate_subdaily_production(control *, fluxes *, met *m, params *,
                                     state *, int, double);

void   calc_autotrophic_respiration(control *, fluxes *, met *, params *,
                                    state *);
double lloyd_and_taylor(double);

/* N stuff */
int    nitrogen_allocation(control *c, fluxes *, params *, state *, double,
                           double, double, double, double, double, int);
double calculate_growth_stress_limitation(params *, state *);
double calculate_nuptake(control *, params *, state *);

double nitrogen_retrans(control *, fluxes *, params *, state *,
                        double, double, int);
void   calculate_ncwood_ratios(control *c, params *, state *, double, double *,
                              double *, double *, double *);

/* Priming/Exudation stuff */
void   calc_root_exudation(control *c, fluxes *, params *p, state *);

/* hydraulics */
void   initialise_roots(fluxes *, params *, state *);
void   update_roots(control *, params *, state *);
double calc_root_dist(double, double, double, double, double, double);



#endif /* PLANT_GROWTH */
