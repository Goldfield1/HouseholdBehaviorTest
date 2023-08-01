typedef struct sol_struct
{
 double* c;
 double* h;
 double* V;
 double* c_egm;
 double* h_egm;
 double* V_egm;
} sol_struct;

double* get_double_p_sol_struct(sol_struct* x, char* name){

 if( strcmp(name,"c") == 0 ){ return x->c; }
 else if( strcmp(name,"h") == 0 ){ return x->h; }
 else if( strcmp(name,"V") == 0 ){ return x->V; }
 else if( strcmp(name,"c_egm") == 0 ){ return x->c_egm; }
 else if( strcmp(name,"h_egm") == 0 ){ return x->h_egm; }
 else if( strcmp(name,"V_egm") == 0 ){ return x->V_egm; }
 else {return NULL;}

}


