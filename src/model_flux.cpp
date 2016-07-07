#include <vector>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include "arc.h"
#include "t.h"

#include "model_flux.h"

using namespace std;
ILOSTLBEGIN

string modelFlux(int l, int N[],int  ord[], bool extra, vector<vector<Arc>> O, int n, int o, int initialSolution[][3]){

    IloEnv env;

    stringstream logfile;
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
        model.add(obj <= l);
        model.add(obj >= (l*2/3));

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
            expr.end();
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
                    expr.end();
                }

            }
        }

        //strips

        int stripInicial = -1;
        if(ord[0] == 0){
            stripInicial = 0;
            while((stripInicial < (n-1)) && N[stripInicial + 1] - N[stripInicial] == 1){
                stripInicial++;
            }

        }

        int stripFinal = n;
        if(ord[n-1] == (n-1)){
            stripFinal = n-1;
            while((stripFinal > 0) && N[stripFinal] - N[stripFinal - 1] == 1){
                stripFinal--;
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
                        bool validOperation = true;
                        for(int aux = 0; aux <= stripInicial; aux++){
                            Arc arcAux(aux,aux);
                            if(!contains(arcAux, O[y])){
                                validOperation = false;
                                break;
                            }
                        }

                        for(int aux = n-1; aux >= stripFinal; aux--){
                            Arc arcAux(aux,aux);
                            if(!contains(arcAux, O[y])){
                                validOperation = false;
                                break;
                            }
                        }
                        //validOperation = true;
                        if(extra == false){

                            validOperation = true;
                        }
                        if(validOperation){
                            Arc arc(a,b);

                            if(contains(arc, O[y])){
                                flowCapacityTerm2 += z[k][y];
                            }
                        }
                    }

                    model.add(flowCapacityTerm1 <= flowCapacityTerm2);

                    flowCapacityTerm1.end();
                    flowCapacityTerm2.end();
                }
            }
        }

        // constraint to u0 be always in the end
        for(int k = 0; k < l-1; k++){
            model.add(z[k][0] <= z[k+1][0]);
        }

        //extra
        if(extra){
            cout << "Restrições extras" << endl;
            for (int k = 0; k < l-1; k++){

                //cout << "k = " << k << endl;
                for (int i = 0; i < (n-1); i++){

                    //cout << "a2 = " << a2 << endl;
                    for(int b = 0; b < (n-1); b++){

                        //cout << "b = " << b << endl;
                            //for(int j = 1; (k+j) <= (l-1); j++){
                                //cout << "b2 = " << b2 << endl;
                                //cout << ord[i] << endl;
                                //cout << ord[i+1] << endl;
                                //int j = (l-1) - k;
                                IloExpr expr1(env);
                                IloExpr expr2(env);
                                for(int a = 0; a < n; a ++){
                                    expr1 += x[k][ord[i]][a][b];
                                    expr2 += x[k][ord[i+1]][a][b+1];
                                }

                                for(int b2 = 0; b2 < (n-1); b2++){
                                    model.add(expr1 + expr2 + x[k+1][ord[i]][b][b2] - x[k+1][ord[i+1]][b+1][b2+1] <= 2);
                                }

                                // expr1 + expr2 + z_invalido <= 2

                                expr1.end();
                                expr2.end();

                                //model.add(expr1 + expr2 + x[k+j][ord[i]][b][b2] <= x[k+j][ord[i+1]][b+1][b2+1] + 2);

                                //model.add(expr1 + expr2 + x[k+1][ord[i]][b][b2] <= x[k+1][ord[i+1]][b+1][b2+1] + 2);
                            //}

                    }
                }
            }

            if(ord[0] == 0){

                for(int i = 0; i <= stripInicial; i++){
                    for(int k = 0; k < l; k++){
                       model.add(x[k][i][i][i] == 1);
                    }
                }
            }

            if(ord[n-1] == (n-1)){

                for(int i = n-1; i >= stripFinal; i--){
                   for(int k = 0; k < l; k++){
                       model.add(x[k][i][i][i] == 1);
                    }
                }
            }

            for(int k = 0; k < (l-1); k++){
                for(int a = 0; a < n; a++){
                    model.add(x[k][ord[0]][a][0] <= x[k+1][ord[0]][0][0]);
                }
            }

            for(int k = 0; k < (l-1); k++){
                for(int a = 0; a < n; a++){
                    model.add(x[k][ord[n-1]][a][n-1] <= x[k+1][ord[n-1]][n-1][n-1]);
                }
            }

            for(int i = 0; i < (n-1); i++){

                if(N[i + 1] - N[i] == 1){

                    for(int b = 0; b < (n-1); b++){

                        model.add(x[0][i][i][b] <= x[0][i+1][i+1][b+1]);

                    }

                }

            }

        }

        //Solving the problem

        //env.setOut(logfile);
        IloCplex cplex(model);
        //cplex.setOut(logfile);
        //cplex.setWarning(logfile);
        //cplex.setError(logfile);

        //Adding initial solution
        IloNumVarArray startVar(env);
        IloNumArray startVal(env);
        int permutation[n];
        int nextPermutation[n];
        for(int i = 0; i < n; i++){
            nextPermutation[i] = i;
        }
        for(int k = 0; k < l; k++){
            vector<Arc> transposicao = gerarTransposicao(initialSolution[k][0],initialSolution[k][1],initialSolution[k][2],n);
            for(int i = 0; i < n; i++){
                  permutation[i] = nextPermutation[i];
            }
            for(int i = 0; i < n; i++){
                for(int a = 0; a < n; a++){
                    for(int b = 0; b < n; b++){

                        Arc arc(a,b);

                        if(contains(arc, transposicao) && i == permutation[a]){
                            nextPermutation[b] = permutation[a];
                            //cout << a << b << permutation[a] << endl;
                            startVar.add(x[k][i][a][b]);
                            startVal.add(1);
                        }else{
                            startVar.add(x[k][i][a][b]);
                            startVal.add(0);
                        }

                    }
                }
            }

        }


        for(int k = 0; k < l; k++){
            for(int y = 0; y < o; y++){
                int finalIndex = 0;
                int index = 0;
                if(initialSolution[k][0] == 0 && initialSolution[k][1] == 0 && initialSolution[k][2] == 0){
                    finalIndex = 0;
                }else{
                    for(int i = 0; i <= (n-2); i++){
                        for(int j = i + 1; j <= (n-1); j++){
                           for(int k2 = j + 1; k2 <= n; k2++){
                               index++;
                               if(i == initialSolution[k][0] && j == initialSolution[k][1] && k2 == initialSolution[k][2]){
                                    finalIndex = index;
                               }
                           }
                        }
                    }
                }

                if(y == finalIndex){
                    startVar.add(z[k][y]);
                    startVal.add(1);
                }else{
                        startVar.add(z[k][y]);
                    startVal.add(0);
                }
            }
        }

        cplex.addMIPStart(startVar, startVal);
        startVal.end();
        startVar.end();

        //cplex.extract(model);

        //timeout
        cplex.setParam(IloCplex::Param::TimeLimit,7200);
        //IloCplex::Param::MIP::Strategy::HeuristicFreq
        //cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq,1);
        //IloCplex::Param::MIP::Display

        //cplex.setParam(IloCplex::Param::MIP::Display,4);
        //cplex.setParam(IloCplex::RootAlg,6);
        //cplex.setParam(IloCplex::NodeAlg,6);
        cplex.setParam(IloCplex::Param::MIP::Strategy::VariableSelect, 3);
        //cplex.setParam(IloCplex::MIPEmphasis, 2);
        //IloCplex::Param::MIP::Limits::RepairTries
        //cplex.setParam(IloCplex::Param::MIP::Limits::RepairTries, 5);
        cplex.setParam(IloCplex::Param::MIP::Tolerances::UpperCutoff, l);
        cplex.setParam(IloCplex::Param::MIP::Tolerances::LowerCutoff, l*2/3);

        //cplex.exportModel("/home/sandro/model1.lp");
        if (cplex.solve()) {
            cout << "Optimal value: " << cplex.getObjValue() << endl;
            cout << "l = " << l << endl;
            cout << cplex.getTime() << endl;

            /*for(int k = 0; k < l; k++){

                for(int i = 0; i < n; i++){
                    cout << "k = " << k << ", i = " << i << endl;
                    for(int a = 0; a < n; a++){

                        for(int b = 0; b < n; b++){
                            cout << cplex.getValue(x[k][i][a][b]) << " ";
                        }

                    }


                    cout << endl;
                }

            }*/

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
    const string str = logfile.str();
    //cout << "Saída " <<  str << endl;

    return str;
}
