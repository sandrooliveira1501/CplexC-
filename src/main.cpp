#include <vector>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include "arc.h"
#include "t.h"
#include "model_flux.h"
#include "model.h"
#include "solver.h"
#include "lagrange.h"
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

    exec("trans", "def", N, n);
	
    //llbp(l,N,O,n,O.size());
	*/

}
