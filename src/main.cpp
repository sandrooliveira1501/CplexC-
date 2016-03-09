#include <vector>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include "arc.h"
#include "t.h"
#include "model_flux.h"
#include "model.h"
#include "solver.h"
#include <fstream>

using namespace std;

int main(int argc, const char *argv[]){
    if(argc < 2){
        cout << "Especifique caminho do arquivo de entrada" << endl;
        return 1;
    }

    const char* filename = argv[1];

    ifstream infile;
    infile.open(filename);

    std::string line;


    while (getline(infile, line)){
        int n = 0,tmp;

        int *N;

        istringstream tmp_iss(line);

        cout << line << endl;

        while (tmp_iss >> tmp) {
            n++;
        }

        N = new int[n];
        tmp = 0;

        istringstream iss(line);
        while(iss >> N[tmp]){
            tmp++;
        }

        getline(infile, line);

        istringstream iss_l(line);
        int l = -1;
        if(!(iss_l >> l)){
            cout << "Entrada inválida" << endl;
            return 1;
        }

        vector<vector<Arc>> O = gerarTransposicoes(n);


        cout << "Multicommodity flow model" << endl;
        modelFlux(l,N,O,n,O.size());


        cout << "Perferct Matching model" << endl;
        model(l,N,O,n,O.size());


        for(int i = 0; i < n; i++){
            N[i] = N[i] + 1;
        }

        cout << "Zanoni and Cid model" << endl;
        exec("trans", "def", N, n);

    }

    infile.close();

   /* int n = 5;
    int N[] = {4,3,2,1,0};
    int l = 6;

    vector<vector<Arc>> O = gerarTransposicoes(n);
    cout << O.size() << endl;

    //modelFlux(l,N,O,n,O.size());

    //model(l,N,O,n,O.size());

    for(int i = 0; i < n; i++){
        N[i] = N[i] + 1;
    }

    exec("trans", "def", N, n);*/

}

/*
int main(){

    IloEnv env;
    try{

        IloModel model(env);
        IloNumVarArray vars(env);
        IloRangeArray c(env);
        IloExpr expr1(env);
        IloExpr expr2(env);

        vars.add(IloNumVar(env, 0.0,40.0));
        vars.add(IloNumVar(env));
        vars.add(IloNumVar(env));

        expr1 -= vars[0];
        expr1 += vars[1];
        expr1 += vars[2];


        c.add(-vars[0] + vars[1] + vars[2] <= 20);
        c.add(vars[0] - 3 * vars[1] + vars[2] <= 30);

        model.add(c);
        model.add(IloMaximize(env, vars[0] + 2*vars[1] + 3*vars[2]));

        IloCplex cplex(model);

        cplex.solve();
        cout << " Max =" << cplex . getObjValue () << endl ;


    }catch(IloException& e){
        cerr << "Concert exception caught: " << e << endl;
    }catch(...){
        cerr << "Unknow exception caught" << endl;
    }

    env.end();

    return 0;
}*/

/*
void model(int l, int N[], vector<vector<Arc>> O, int n, int o){

    IloEnv env;

    try{

        IloModel model(env);

        // var z_{k,mi}
        IloArray<IloBoolVarArray> z(env, l);

        for(int i = 0; i < l; i++){

            z[i] = IloBoolVarArray (env, o);

        }

        // var x_{kiab}
        IloArray<IloArray<IloArray<IloIntVarArray>>> x(env, l);
        for (int i = 0; i < l; i++) {
            x[i] = IloArray<IloArray<IloIntVarArray> > (env, n);
            for (int j = 0; j < n; j++) {
                x[i][j] = IloArray<IloIntVarArray> (env, n);
                for (int k = 0; k < n; k++) {
                    x[i][j][k] = IloIntVarArray (env,n,0,IloIntMax);
                }
            }
        }

        //objective function

        IloExpr obj(env);

        for(int k = 0; k < l;k++){
            for(int i = 1; i < o; i++){
                obj += z[k][i];
            }
        }

        model.add(IloMinimize(env, obj));

        //one operator by level
        for(int k = 0; k < l; k++){
            model.add(IloSum(z[k]) == 1);
        }


        // outgoing
        for(int i = 0; i < n; i++){
           model.add(IloSum(x[0][i][i]) == 1);
        }

        // ingoing
        for(int i = 0; i < n; i++){
            IloExpr expr(env);
            for(int j = 0; j < n; j++){
                expr += x[l-1][i][j][N[i]];
            }
            model.add(expr == 1);
        }

        //conservation
        for(int i = 0; i < n; i++){
            for(int k = 1; k < l; k++){
                for(int a = 0; a < n; a++){
                    IloExpr expr(env);

                    for(int j = 0; j < n; j++){
                        expr +=  x[k][i][a][j];
                        expr -=  x[k-1][i][j][a];
                    }

                    model.add(expr == 0);
                }

            }
        }

        // just the arcs that are in mi can be used


            for(int a = 0; a < n; a++){
                for(int b = 0; b < n; b++){
                    for(int k  = 0; k < l; k++){
                        IloExpr flowCapacityTerm1(env);
                        IloExpr flowCapacityTerm2(env);

                        for(int i = 0; i < n; i++){
                            flowCapacityTerm1 += x[k][i][a][b];
                        }

                        for(int y = 0; y < o; y++){
                            Arc arc(a,b);
                            if(contains(arc, O[y])){
                                flowCapacityTerm2 += z[k][y];
                            }
                        }

                        model.add(flowCapacityTerm1 <= flowCapacityTerm2);
                }
            }
        }


        // constraint to u0 be always in the end
        for(int k = 0; k < l-1; k++){
            model.add(z[k][0] <= z[k+1][0]);
        }

        //Solving the problem
        IloCplex cplex(model);
        cplex.extract(model);
        if (cplex.solve()) {
            for (int k = 0; k < l; k++){
                for(int i = 0; i < n; i++){
                    cout << " k " << k << " i " << i << endl;
                    for(int a = 0; a < n; a++){

                        for(int b = 0; b < n; b++){
                            if(cplex.getValue(x[k][i][a][b]) == 1){
                                cout << "a " << a << " b " << b << " ";
                            }
                        }

                    }
                    cout << endl;

                }

            }
            cplex.out() << "Optimal value: " << cplex.getObjValue() << endl;
        }


        obj.end();


    }catch(IloException& e){
        cerr << "Concert exception caught: " << e << endl;
    }catch(...){
        cerr << "Unknow exception caught" << endl;
    }
        env.end();

}*/

