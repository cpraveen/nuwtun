#include "IpStdCInterface.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#define NPMAX 1000
#define NDMAX 50

// Database to store information on computed solutions
typedef struct db {
   Index np;
   Index nd;
   Number x[NPMAX][NDMAX];
   char dir[NPMAX][48];
   } DB;

// Structure to hold input parameters read from file
typedef struct problem {
   Index nd;          // number of design variables
   Index nc;          // number of constraints
   char obj[48];      // objective function
   char con[100][48]; // constraints
   Number *xl, *xu;   // lower and upper bounds
   } Problem;
   
typedef struct result {
   Number cl;
   Number cd;
   Number apgrad;
   Number area;
   } Result;
   
/* Function Declarations */
Bool eval_f(Index n, Number* x, Bool new_x,
            Number* obj_value, UserDataPtr user_data);

Bool eval_grad_f(Index n, Number* x, Bool new_x,
                 Number* grad_f, UserDataPtr user_data);

Bool eval_g(Index n, Number* x, Bool new_x,
            Index m, Number* g, UserDataPtr user_data);

Bool eval_jac_g(Index n, Number *x, Bool new_x,
                Index m, Index nele_jac,
                Index *iRow, Index *jCol, Number *values,
                UserDataPtr user_data);

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nele_hess, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data);

void init();
void finddir(Number *x, char *dir);
void add2db(Number *x, char *dir);
void read_result(char *dir, Result *res);
void read_grad(char *dir, Index n, Number *grad_f);
void solve_flo(Index n, Number *x, Result *res);

/* Global variables */
DB optdb;
Number clref, cdref, apgradref, arearef;
Problem param;
Result res0;

char mkdir[48];
char hicks[48];
char deform[48];
char deform_adj[48];
char nuwtun_flo[48];
char nuwtun_adj[48];
char nuwtun_area[48];
char nuwtun_areax[48];
char refval[48];

/* Main Program */
int main()
{
  Index n=-1;                          /* number of variables */
  Index m=-1;                          /* number of constraints */
  Number* x_L = NULL;                  /* lower bounds on x */
  Number* x_U = NULL;                  /* upper bounds on x */
  Number* g_L = NULL;                  /* lower bounds on g */
  Number* g_U = NULL;                  /* upper bounds on g */
  IpoptProblem nlp = NULL;             /* IpoptProblem */
  enum ApplicationReturnStatus status; /* Solve return code */
  Number* x = NULL;                    /* starting point and solution vector */
  Number* mult_x_L = NULL;             /* lower bound multipliers
                                          at the solution */
  Number* mult_x_U = NULL;             /* upper bound multipliers
                                          at the solution */
  Number obj;                          /* objective value */
  Index i;                             /* generic counter */
  char dir[48];

  init(); // Read parameters from ipopt.in file

  /* set the number of variables and allocate space for the bounds */
  n=param.nd;
  x_L = (Number*)malloc(sizeof(Number)*n);
  x_U = (Number*)malloc(sizeof(Number)*n);
  /* set the values for the variable bounds */
  for (i=0; i<n; i++) {
    x_L[i] = param.xl[i];
    x_U[i] = param.xu[i];
  }

  /* set the number of constraints and allocate space for the bounds */
  m = param.nc;
  g_L = (Number*)malloc(sizeof(Number)*m);
  g_U = (Number*)malloc(sizeof(Number)*m);
  /* set the values of the constraint bounds */
  for(i=0; i<m; i++){
     g_L[i] = 0.0;
     g_U[i] = 0.0;
  }

  /* initialize database */
  optdb.np = 0;
  optdb.nd = n;

  /* Number of nonzeros in the Jacobian of the constraints */
  Index nele_jac = n*m;
  /* Number of nonzeros in the Hessian of the Lagrangian (lower or
     upper triangual part only) */
  Index nele_hess = n*(n+1)/2;
  /* indexing style for matrices */
  Index index_style = 0; /* C-style; start counting of rows and column
                            indices at 0 */

  /* create the IpoptProblem */
  nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
                           index_style, &eval_f, &eval_g, &eval_grad_f,
                           &eval_jac_g, &eval_h);

  /* We can free the memory now - the values for the bounds have been
     copied internally in CreateIpoptProblem */
  free(x_L);
  free(x_U);
  free(g_L);
  free(g_U);

  /* Set some options.  Note the following ones are only examples,
     they might not be suitable for your problem. */
  AddIpoptNumOption(nlp, "tol", 1e-2);
  AddIpoptNumOption(nlp, "constr_viol_tol", 1e-2);
  AddIpoptNumOption(nlp, "acceptable_tol", 1e-2);
  AddIpoptNumOption(nlp, "dual_inf_tol", 1e-2);
  AddIpoptNumOption(nlp, "compl_inf_tol", 1e-2);
  AddIpoptStrOption(nlp, "mu_strategy", "monotone");
  AddIpoptStrOption(nlp, "output_file", "ipopt.out");
  AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
  AddIpoptIntOption(nlp, "limited_memory_max_history", 6);
  AddIpoptStrOption(nlp, "print_user_options", "yes");
  AddIpoptIntOption(nlp, "print_level", 5);
  AddIpoptIntOption(nlp, "max_iter", 30);

  /* allocate space for the initial point and set the values */
  x = (Number*)malloc(sizeof(Number)*n);
  for(i=0; i<n; i++) x[i] = 0.0;

  /* allocate space to store the bound multipliers at the solution */
  mult_x_L = (Number*)malloc(sizeof(Number)*n);
  mult_x_U = (Number*)malloc(sizeof(Number)*n);

  /* delete old workdirs */
  system("rm -rf workdir.*");

  /* solve the problem */
  status = IpoptSolve(nlp, x, NULL, &obj, NULL, mult_x_L, mult_x_U, NULL);

  finddir(x, dir);
  printf("Best solution is in directory %s\n", dir);

  //if (status == Solve_Succeeded) {
  if (1) {
    printf("\n\nSolution of the primal variables, x\n");
    for (i=0; i<n; i++) {
      printf("x[%d] = %e\n", i, x[i]);
    }

    printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
    for (i=0; i<n; i++) {
      printf("z_L[%d] = %e\n", i, mult_x_L[i]);
    }
    for (i=0; i<n; i++) {
      printf("z_U[%d] = %e\n", i, mult_x_U[i]);
    }

    printf("\n\nObjective value\n");
    printf("f(x*) = %e\n", obj);
  }

  /* free allocated memory */
  FreeIpoptProblem(nlp);
  free(x);
  free(mult_x_L);
  free(mult_x_U);

  return 0;
}

//=============================================================================
// Read parameters from file ipopt.in
//=============================================================================
void init(){
   FILE *fpt;
   int i;
   
   fpt = fopen("ipopt.in", "r");
   fscanf(fpt,"%s", param.obj);
   fscanf(fpt,"%d", &param.nc);
   for(i=0; i<param.nc; i++) fscanf(fpt,"%s", param.con[i]);
   fscanf(fpt,"%d", &param.nd);
   param.xl = (Number*)malloc(sizeof(Number)*param.nd);
   param.xu = (Number*)malloc(sizeof(Number)*param.nd);
   for(i=0; i<param.nd; i++)
      fscanf(fpt,"%lf%lf", &param.xl[i], &param.xu[i]);
   fclose(fpt);
   
   printf("Objective function    = %s\n", param.obj);
   printf("Number of constraints = %d\n", param.nc);
   for(i=0; i<param.nc; i++) 
      printf("  Constraint[%d]       = %s\n", i+1, param.con[i]);
   printf("Number of design vars = %d\n", param.nd);
   for(i=0; i<param.nd; i++)
      printf("  Variable[%2d] = %12.4e  %12.4e\n", 
             i+1, param.xl[i], param.xu[i]);

}
//=============================================================================
/* Objective function */
//=============================================================================
Bool eval_f(Index n, Number* x, Bool new_x,
            Number* obj_value, UserDataPtr user_data)
{
  Result res;
    
  solve_flo(n, x, &res);

  if(strcmp(param.obj,"CD")==0){
     *obj_value = res.cd/res0.cd;
  }else
  if(strcmp(param.obj,"APGRAD")==0){
     *obj_value = res.apgrad/res0.apgrad;
  }else{
     printf("eval_f: Unknown objective function !!!\n");
     exit(0);
  }

  return TRUE;
}

//=============================================================================
/* Gradient of objective function */
//=============================================================================
Bool eval_grad_f(Index n, Number* x, Bool new_x,
                 Number* grad_f, UserDataPtr user_data)
{
  char dir[48];
  Result res;

  finddir(x, dir);
  if(!strcmp(dir,"")){
     solve_flo(n, x, &res);
     finddir(x, dir);
  }
  sprintf(refval,"cp workdir.1/fort.19 %s/refval.dat", dir);
  system(refval);
  sprintf(nuwtun_adj,"cd %s && nuwtun_adj %s < flo.in > adj_cd.log", 
          dir, param.obj);
  system(nuwtun_adj);
  sprintf(deform_adj,"cd %s && deform_adj > def_adj.log", dir);
  system(deform_adj);
  read_grad(dir, n, grad_f);

  return TRUE;
}

//=============================================================================
/* Constraint functions */
//=============================================================================
Bool eval_g(Index n, Number* x, Bool new_x,
            Index m, Number* g, UserDataPtr user_data)
{
  Result res;
  int i;
  
  solve_flo(n, x, &res);

  for(i=0; i<m; i++){
     if(strcmp(param.con[i],"CL")==0){
        g[i] = res.cl/res0.cl - 1.0;
     }else if(strcmp(param.con[i],"AREA")==0){
        g[i] = res.area/res0.area - 1.0;
     }else{
        printf("eval_g: Unknown constraint function !!!\n");
        exit(0);
     }
  }

  return TRUE;
}

//=============================================================================
/* Jacobian of constraints */
//=============================================================================
Bool eval_jac_g(Index n, Number *x, Bool new_x,
                Index m, Index nele_jac,
                Index *iRow, Index *jCol, Number *values,
                UserDataPtr user_data)
{
  char dir[48];
  Index i, c=0;
  Result res;
  Number *grad_g;

  grad_g = (Number*)malloc(sizeof(Number)*n);

  if (values == NULL) {
    /* return the structure of the jacobian */
    /* this particular jacobian is dense */
     int i, j, c = 0;
     for(i=0; i<m; i++)
        for(j=0; j<n; j++){
           iRow[c] = i;
           jCol[c] = j;
           ++c;
        }
  }
  else {
    /* return the values of the jacobian of the constraints */

      finddir(x, dir);
      if(!strcmp(dir,"")){
         solve_flo(n, x, &res);
         finddir(x, dir);
      }
      sprintf(refval,"cp workdir.1/fort.19 %s/refval.dat", dir);
      system(refval);
   
      for(i=0; i<m; i++){
         
         if(strcmp(param.con[i],"CL")==0){
            sprintf(nuwtun_adj,"cd %s && nuwtun_adj %s < flo.in > adj_cl.log", 
                    dir, param.con[i]);
            system(nuwtun_adj);
         }
         else if(strcmp(param.con[i],"AREA")==0){
            sprintf(nuwtun_areax,"cd %s && area_x > adj_area.log", 
                    dir);
            system(nuwtun_areax);
         }else{
            printf("eval_jac_g: Unknown constraint function !!!\n");
            exit(0);
         }
         
         sprintf(deform_adj,"cd %s && deform_adj > def_adj.log", dir);
         system(deform_adj);
         read_grad(dir, n, grad_g);
      
         /* For area, scale gradient by reference area */
         if(strcmp(param.con[i],"AREA")==0)
            for(i=0; i<n; i++)
               grad_g[i] = grad_g[i]/res0.area;

         for(i=0; i<n; i++)
          values[c++] = grad_g[i];
      }

  }

  free(grad_g);
  return TRUE;
}

//=============================================================================
/* Hessian of Lagrangian function. Not required since we use quasi-newton */
//=============================================================================
Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nele_hess, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data)
{
   return FALSE;
}

//=============================================================================
// Find if x exists in database. If yes then return its directory
//=============================================================================
void finddir(Number *x, char *dir){
   Index i, j;
   int stat;
   sprintf(dir, "%s", ""); // set to empty string
   for(i=0; i<optdb.np; i++){
      stat = 1;
      for(j=0; j<optdb.nd; j++)
         if(optdb.x[i][j] != x[j]) stat = 0;
      if(stat==1){
         sprintf(dir, "%s", optdb.dir[i]);
         return;
      }
   }
}

//=============================================================================
// Add point x to database and set workdir
//=============================================================================
void add2db(Number *x, char *dir){
   Index i, j, dirnum;

   if(optdb.np == NPMAX || optdb.nd > NDMAX){
      printf("add2db: database size is small\n");
      printf("NPMAX = %d, required is >  %d\n", NPMAX, optdb.np);
      printf("NDMAX = %d, required is >= %d\n", NDMAX, optdb.nd);
      exit(0);
   }

   i = optdb.np;
   for(j=0; j<optdb.nd; j++) optdb.x[i][j] = x[j];

   dirnum = i + 1;
   sprintf(dir,"%s%d", "workdir.", dirnum);
   //printf("i=%d, dir=%s\n",i,dir);
   sprintf(optdb.dir[i],"%s%d", "workdir.", dirnum);
   optdb.np += 1;
}

//=============================================================================
// Read cl and cd from file
//=============================================================================
void read_result(char *dir, Result *res){
   FILE *fpt;
   char str[48], file[48];
   Index idum;
   Number fdum;

   sprintf(file,"%s/fort.19", dir);
   fpt = fopen(file, "r");
   if(fpt == NULL){
      printf("read_result: could not open %s\n", file);
      exit(0);
   }
   fscanf(fpt,"%s%d", str, &idum);
   fscanf(fpt, "%s%lf%lf%lf%lf%lf", str, &fdum, &fdum, &fdum, &fdum, &fdum);
   fscanf(fpt, "%s%lf", str, &(res->cl));
   fscanf(fpt, "%s%lf", str, &(res->cd));
   fscanf(fpt, "%s%lf", str, &(res->apgrad));
   fclose(fpt);

   /* Read area.out if it exists */
   sprintf(file,"%s/area.out", dir);
   fpt = fopen(file, "r");
   if(fpt == NULL){
      res->area = 0.0;
   }else{
      fscanf(fpt, "%lf", &(res->area));
   }
   fclose(fpt);

}

//=============================================================================
/* Run flow solver */
//=============================================================================
void solve_flo(Index n, Number *x, Result *res)
{
  char dir[48];
  FILE *fpt;
  Index i;
  
  finddir(x, dir);
  if(!strcmp(dir,"")){
     add2db(x,dir);
     sprintf(mkdir,"cp -r templatedir %s", dir);
     system(mkdir);
     sprintf(hicks,"%s/hicks.in", dir);
     fpt = fopen(hicks, "w");
     for(i=0; i<n; i++) fprintf(fpt,"%20.10e\n", x[i]);
     fclose(fpt);
     sprintf(deform,"cd %s && deform > def.log", dir);
     system(deform);
     sprintf(nuwtun_flo,"cd %s && nuwtun_flo < flo.in > flo.log", dir);
     system(nuwtun_flo);
     for(i=0; i<param.nc; i++) //If area is a constraint, then compute it
        if(strcmp(param.con[i],"AREA")==0){
           sprintf(nuwtun_area,"cd %s && area > area.log", dir);
           system(nuwtun_area);
        }
  }
  read_result(dir, res);
  if(!strcmp(dir,"workdir.1")){
     res0.cl     = res->cl;
     res0.cd     = res->cd;
     res0.apgrad = res->apgrad;
     res0.area   = res->area;
  }

}

//=============================================================================
// Read gradient from gradient.dat
//=============================================================================
void read_grad(char *dir, Index n, Number *grad_f){
   char gradfile[48];
   FILE *fpt;
   Index i, j;

   sprintf(gradfile, "%s/gradient.dat", dir);
   fpt = fopen(gradfile, "r");
   if(fpt == NULL){
      printf("read_grad: could not open %s\n", gradfile);
      exit(0);
   }
   for(i=0; i<n; i++) fscanf(fpt,"%d%lf", &j, &grad_f[i]);
   fclose(fpt);
}
