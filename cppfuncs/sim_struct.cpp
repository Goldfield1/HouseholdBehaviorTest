typedef struct sim_struct
{
 double* c;
 double* h;
 double* a;
 double* k;
 double* a_init;
 double* k_init;
} sim_struct;

double* get_double_p_sim_struct(sim_struct* x, char* name){

 if( strcmp(name,"c") == 0 ){ return x->c; }
 else if( strcmp(name,"h") == 0 ){ return x->h; }
 else if( strcmp(name,"a") == 0 ){ return x->a; }
 else if( strcmp(name,"k") == 0 ){ return x->k; }
 else if( strcmp(name,"a_init") == 0 ){ return x->a_init; }
 else if( strcmp(name,"k_init") == 0 ){ return x->k_init; }
 else {return NULL;}

}


