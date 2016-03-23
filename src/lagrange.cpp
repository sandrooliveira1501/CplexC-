#include "lagrange.h"

using namespace std;
ILOSTLBEGIN

Lagrange::Lagrange(int l, int *N, vector<vector<Arc>> O, int n, int o, int LB){

   for(int i = 0; i < l; i++){
        alpha[i] = -1;
        for(int a = 0; a < n; a++){
            for(int b = 0; b < n; b++){
                beta[i][a][b] = -1;
            }
        }
    }

    this->LB = LB;
    this->ZUB = l;
    this->ZLB = 0;
    this->l = l;
    this->N = N;
    this->O = O;
    this->n = n;
    this->o = o;
}

void Lagrange::execute(){

    IloEnv env;
    float pi = 2;
    float auxZLB = 0;
    try{
        int n = 0;
        while(pi > 0.005){
            IloModel model = prepareModel(env);

            env.setOut(env.getNullStream());

            IloCplex cplex(model);

            cplex.setOut(env.getNullStream());
            cplex.setWarning(env.getNullStream());
            cplex.setError(env.getNullStream());

            cplex.extract(model);
            cplex.exportModel("/home/alexsandro/llp.lp");

            if (cplex.solve()) {

                cout << "Optimal value: " << cplex.getObjValue() << endl;
                cout << cplex.getTime() << endl;

                cplex.getObjValue();
                this->ZLB = (cplex.getObjValue() > this->ZLB)? cplex.getObjValue(): this->ZLB;

                //transformar ZLB em uma solução viável

                //calcular G e T

                int G[this->l];
                int G2 [this->l][this->n][this->n];
                int aux = 0;
                int aux2 = 0;

                for(int k = 0; k < this->l; k++){

                    G[k] = 1;

                    for(int o = 0; o < this->o; o++){
                        G[k] -= cplex.getValue(this->z[k][o]);
                    }

                    aux += (G[k] * G[k]);
                }

                for(int k = 0; k < this->l; k++){

                    for(int a = 0; a < this->n; a++){

                        for(int b = 0; b < this->n; b++){

                            G2[k][a][b] = 0;

                            for(int i = 0; i < this->n; i++){

                                G2[k][a][b] += cplex.getValue(x[k][i][a][b]);

                            }

                            for(int o = 0; o < this->o; o++){

                                Arc arc(a,b);
                                if(contains(arc, this->O[o])){
                                    G2[k][a][b] -= cplex.getValue(z[k][o]);
                                }

                            }

                            aux2 += (G2[k][a][b] * G2[k][a][b]);

                        }

                    }

                }


                float T = pi * (1.05*this->ZUB - this->ZLB)/aux;

                float T2 = pi * (1.05*this->ZUB - this->ZLB)/aux2;

                //recalcular alpha e beta

                for(int k = 0; k < this->l; k++){
                    this->alpha[k] = this->alpha[k] + T*G[k];
                }

                for(int k = 0; k < this->l; k++){

                    for(int a = 0; a < this->n; a++){

                        for(int b = 0; b < this->n; b++){

                            float aux = this->beta[k][a][b] + T2 * G2[k][a][b];

                            this->beta[k][a][b] = (aux > 0) ? aux : 0;

                        }

                    }

                }

                if(this->ZLB <= auxZLB){
                    n+=1;
                }else{
                    auxZLB = this->ZLB;
                    n = 0;
                }

                if(n == 30){
                    n = 0;
                    pi /= 2;
                }
            }else{
                cout << "timeout" << endl;
            }

            cplex.end();
            model.end();

        }
    }catch(IloException& e){
        cerr << "Concert exception caught: " << e << endl;
    }catch(...){
        cerr << "Unknow exception caught" << endl;
    }

    env.end();
}



IloModel Lagrange::prepareModel(IloEnv env){

    IloModel model(env);

    z = IloArray<IloBoolVarArray>(env, l);

    for(int i = 0; i < l; i++){

        z[i] = IloBoolVarArray (env, o);

    }

    x = IloArray<IloArray<IloArray<IloIntVarArray>>>(env, l);
    for (int i = 0; i < l; i++) {
        x[i] = IloArray<IloArray<IloIntVarArray> > (env, n);
        for (int j = 0; j < n; j++) {
            x[i][j] = IloArray<IloIntVarArray> (env, n);
            for (int k = 0; k < n; k++) {
                x[i][j][k] = IloIntVarArray (env,n,0,INFINITY);
            }
        }
    }

    IloExpr obj(env);

    for(int k = 0; k < l;k++){
        for(int u = 0; u < o; u++){

            float c;

            c = 1 - alpha[k];

            for(int e = 0; e < n; e++){
                int a = O[u][e].a;
                int b = O[u][e].b;

                c -= beta[k][a][b];

            }

            obj+= z[k][u]*c;

        }

        obj += -z[k][0];

        obj += alpha[k];

        for(int i = 0; i < n; i++){

            for(int a = 0; a < n; a++){

                for(int b = 0; b < n; b++){

                    obj += beta[k][a][b]*x[k][i][a][b];

                }

            }

        }

    }

    model.add(IloMinimize(env, obj));
    model.add(obj >= this->LB);

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

    //capacite
    for(int k = 0; k < l; k++){
        for(int a = 0; a < n; a++){
            for(int b = 0; b < n; b++){
                IloExpr expr(env);

                for(int i = 0; i < n; i++){

                    expr += x[k][i][a][b];

                }

                model.add(expr <= 1);

            }
        }
    }

    return model;
}
