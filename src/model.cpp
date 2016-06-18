#include <vector>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include "arc.h"
#include "t.h"
#include "model.h"

using namespace std;
ILOSTLBEGIN

void model(int l, int N[], std::vector<std::vector<Arc>> O, int n, int o){

    IloEnv env;

    try{

        IloModel model(env);

        IloArray<IloArray<IloBoolVarArray>> b(env,n);

        for(int i = 0; i < n; i++){
            b[i] = IloArray<IloBoolVarArray> (env, n);

            for(int j = 0; j < n; j++){
                b[i][j] = IloBoolVarArray(env,l);
            }
        }


        IloArray<IloBoolVarArray> m(env,o);
        for(int t = 0; t < o; t++){
            m[t] = IloBoolVarArray(env,l);
        }

        IloArray<IloNumVarArray> v(env, l+1);
        for(int k = 0; k < l+1;k++){
            v[k] = IloNumVarArray(env,n,0,n-1);
        }

        IloExpr obj(env);
        for(int t = 1; t < o; t++){
            for(int k = 0; k < l; k++){
                obj += m[t][k];
            }
        }

        model.add(IloMinimize(env, obj));

        for(int a = 0; a < n; a++){
            model.add(v[0][a] == N[a]);
        }

        for(int a = 0; a < n; a++){
            model.add(v[l][a] == a);
        }

        for(int k = 1; k < (l+1); k++){
            for(int a = 0; a < n; a++){
                for(int i = 0; i < n; i++){
                    IloIfThen expr(env, b[i][a][k-1] == 1, v[k][a] == v[k-1][i]);
                    //expr += (b[i][a][k-1] * v[k-1][i]);
                    model.add(expr);
                    //problem ending IloIfThen expression
                }
            }
        }


        for(int k = 0; k < l; k++){
            IloExpr expr(env);

            for(int t = 0; t < o; t++){
                expr += m[t][k];
            }

            model.add(expr == 1);
            expr.end();
        }


        for(int k = 0; k < l; k++){

            for(int j = 0; j < n; j++){
                IloExpr expr(env);
                for(int i = 0; i < n; i++){
                    expr += b[i][j][k];
                }
                model.add(expr == 1);
                expr.end();
            }

        }


        for(int k = 0; k < l; k++){

            for(int i = 0; i < n; i++){
                IloExpr expr(env);
                for(int j = 0; j < n; j++){
                    expr += b[i][j][k];
                }
                model.add(expr == 1);
                expr.end();
            }

        }



        for(int k = 0; k < l-1; k++){

            model.add(m[0][k] <= m[0][k+1]);

        }

        for(int k = 0; k < l; k++){
            for(int t = 0; t < o; t++){
                IloExpr expr(env);

                vector<Arc> transposition =  O[t];

                for(int x = 0; x < n; x++){
                    Arc arc = transposition[x];
                    expr += b[arc.a][arc.b][k];
                 }

                IloIfThen ifExpr(env, m[t][k] == 1, expr == n);
                model.add(ifExpr);
                expr.end();
                //problem ending IloIfThen expression
            }

        }

        //solving the problem
        //env.setOut(env.getNullStream());
        IloCplex cplex(model);
        //cplex.setOut(env.getNullStream());
        //cplex.setWarning(env.getNullStream());
        //cplex.setError(env.getNullStream());
        //cplex.extract(model);

        //cplex.exportModel("/home/alexsandro/model2.lp");
        cplex.setParam(IloCplex::Param::TimeLimit,7200);

        if (cplex.solve()) {
            cout << "Optimal value: " << cplex.getObjValue() << endl;
            cout << cplex.getTime() << endl;
        }else{
            cout << "timeout" << endl;
        }

        cplex.end();
        model.end();
        obj.end();


    }catch(IloException& e){
        cerr << "Concert exception caught: " << e << endl;
    }catch(...){
        cerr << "Unknow exception caught" << endl;
    }

    env.end();

}
