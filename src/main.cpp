//
// Created by guillermo on 6/11/18.
//
#include <algorithm>
#include <numeric>

#include "problemcec2014.h"
#include "random.h"
#include "srandom.h"
#include "YYOPBase.h"
#include "AYYOPBase.h"
#include "YYOPMeme1.h"
#include "YYOPMemeDif.h"
#include "YYOPMemeDifStop.h"

using namespace realea;

void getInitRandom(Random *random, DomainRealPtr domain, tChromosomeReal &crom) {
    tReal min, max;

    for (unsigned i = 0; i < crom.size(); ++i) {
        domain->getValues(i, &min, &max, true);
        crom[i] = random->randreal(min, max);
    }

}

int main(int argc, char *argv[]) {
    int dim = 30;
    //int fun = 1;
    unsigned max_iter = 3;

    DomainRealPtr domain;
    tChromosomeReal sol(dim);
    ProblemCEC2014 cec2014(dim);

    int seed=13;
    Random random(new SRandom(seed));

    int rounds = 5;
    cout << dim << " " << rounds << endl;

    vector<int> fs = {1,11,17};
    //for(int fun=1;fun<=20;fun++){
    for(int fun: fs){
        // Get the function fun for dimension dim
        ProblemPtr problem = cec2014.get(fun);
        // Domain is useful for clipping solutions
        domain = problem->getDomain();
        vector<tFitness> fList(rounds);
        for(int i=0; i<rounds; i++){
            // Init the initial solution (for LS)
            getInitRandom(&random, domain, sol);
            // Get the maximum evaluations from the problem
            unsigned max_evals = problem->getMaxEval();
            max_iter = max_evals;
            //YYOPBase alg;
            //YYOPMeme1 alg;
            YYOPMemeDifStop alg;

            alg.setRandom(&random);
            alg.setProblem(problem.get());

            tFitness fitness = problem->eval(sol);

            alg.apply(sol,fitness,max_iter);
            fList[i] = fitness;
        }
        auto minmaxResult = minmax_element(fList.begin(),fList.end());
        tFitness best = *minmaxResult.first;
        tFitness worst = *minmaxResult.second;
        tFitness mean = accumulate(fList.begin(),fList.end(),0.0)/rounds;
        cout << "#" << fun << ": "
             << std::scientific << mean << " "
             << std::scientific << best << " "
             << std::scientific << worst << " "
             << std::scientific << problem.get()->getOptime() << endl;
    }
}