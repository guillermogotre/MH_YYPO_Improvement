//
// Created by guillermo on 6/11/18.
//

#include "YYOPBase.h"

#include "localsearch.h"
#include "cmaeshan.h"
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <float.h>

tFitness realea::YYOPBase::evalSol(tChromosomeReal &orig) {
    tChromosomeReal fitted(orig.size());
    std::transform(orig.begin(),orig.end(),fitted.begin(),
                   [&](tReal v){return v*(maxV-minV)+minV;}
    );
    return m_problem->eval(fitted);
}

realea::tChromosomeReal realea::YYOPBase::adaptSol(tChromosomeReal & s){
    tChromosomeReal a(s.size());
    std::transform(s.begin(),s.end(),a.begin(),
                   [&](tReal v){return (v-minV)/(maxV-minV);}
    );
    return a;
}

tFitness realea::YYOPBase::splitStage(tChromosomeReal & s, tReal d){
    const double SPLIT_PROB = 0.5;
    const unsigned D = s.size();
    const double DWAYFACTOR = sqrt(2);

    vector<tChromosomeReal> S(2*D,s);

    double R = m_random->randreal(0,1);
    //One-way splitting
    if(R > SPLIT_PROB){
        for(unsigned i=0;i<s.size();i++){
            //I
            tReal newR  = m_random->randreal(0,1)*d + s[i];
            if(newR > 1 || newR < 0)
                newR = m_random->randreal(0,1);
            S[i][i] = newR;
            //D+I
            newR  = m_random->randreal(0,1)*d - s[i];
            if(newR > 1 || newR < 0)
                newR = m_random->randreal(0,1);
            S[D+i][i] = newR;
        }
        //D-way splitting
        //TODO not equal
    }else{
        for(unsigned i=0; i<S.size(); i++){
            for(unsigned j=0; j<D; j++){
                tReal newR;
                if(m_random->randint(0,1)==1){
                    newR  = m_random->randreal(0,1)*(d/DWAYFACTOR) + s[j];
                }else{
                    newR  = m_random->randreal(0,1)*(d/DWAYFACTOR) - s[j];
                }
                if(newR > 1 || newR < 0)
                    newR = m_random->randreal(0,1);
                S[i][j] = newR;
            }
        }
    }
    int bestIndex = -1;
    tFitness bestFitness = DBL_MAX;
    for(unsigned i=0; i<S.size(); i++){
        tFitness f = evalSol(S[i]);
        if(f<bestFitness){
            bestFitness = f;
            bestIndex = i;
        }
    }
    assert(bestIndex != -1);
    s = S[bestIndex];
    return bestFitness;
}

unsigned realea::YYOPBase::apply(tChromosomeReal &sol, tFitness &fitness, unsigned itera) {

    const unsigned I_MIN = 2, I_MAX = 6, I = m_random->randint(I_MIN,I_MAX);
    const unsigned MAX_ITER = itera;
    const double MAX_D = 0.75;

    //cout << "I: " << I << endl;
    int iters = itera;
    tReal alpha = 10;
    tReal d1 = 0.5, d2 = 0.5;
    tChromosomeReal p1 = adaptSol(sol), p2(sol.size());

    for(unsigned i=0; i<sol.size(); i++){
        p2[i] = m_random->randreal(0,1);
    }

    vector<tChromosomeReal> archive(2*I);
    vector<tReal> archiveEval(2*I);

    tChromosomeReal fittedSol(sol.size());

    tReal eval1 = evalSol(p1);
    tReal eval2 = evalSol(p2);
    int count1 = 0, count2 = 0;
    iters-=2;

    tReal bestEval = DBL_MAX;

    //Aux functions
    auto swp = [&](){p1.swap(p2);swap(eval1,eval2),swap(d1,d2);swap(count1,count2);};

    int i=0;
    int added = p1.size()*2*I;
    while((iters-added) > 0 || i != 0){ //After archive stage
        //Check fittest
        if(eval2 < eval1){
            swp();
            cout << MAX_ITER-iters << " " << i << " " << d1 << " " << d2 << " SWAP!" << " ["<< count1 <<"," <<  count2<<"] "<< endl;
        }
        if(eval1 < bestEval)
            bestEval = eval1;

        //Add to archive
        unsigned offset = i*2;

        archive[offset] = p1;
        archive[offset+1] = p2;
        archiveEval[offset] = eval1;
        archiveEval[offset+1] = eval2;

        i++;
        count1++;

        //Archive stage
        if(i == I){
            i=0;
            //Select two fittest
            vector<unsigned> indices(archive.size());
            unsigned n = 0;
            generate(indices.begin(),indices.end(),[&]{ return n++;});
            sort(indices.begin(),indices.end(),
                 [&](int i1, int i2){
                     return archiveEval[i1] < archiveEval[i2];
                 });
            p1 = archive[indices[0]];
            eval1 = archiveEval[indices[0]];

            p2 = archive[indices[1]];
            eval2 = archiveEval[indices[1]];

            //Update d
            d1 = d1 - (d1/alpha);
            d2 = min(d2 + (d2/alpha),MAX_D);
        }
            //Splitting stage
        else{
            eval1 = splitStage(p1,d1);
            eval2 = splitStage(p2,d2);

            iters -= 4*p1.size();
        }
    }

    /*
    ILocalSearch *ls;
    ILSParameters *ls_options;

    CMAESHansen *cmaes = new CMAESHansen("cmaesinit.par");
    cmaes->searchRange(0.1);

    ls = cmaes;
    ls->setProblem(m_problem);
    ls->setRandom(m_random);

    ls_options = ls->getInitOptions(sol);
    */

    assert(bestEval == eval1);

    tFitness before,after,diff;

    before = fitness;

    //cout << "Evals: " << MAX_ITER - iters << endl;

    after = eval1;

    diff = before-after;

    //cout <<"Improvement: " <<std::scientific <<before <<" -> " <<std::scientific <<after;
    //cout <<" [" <<std::scientific <<diff <<"]" << " , Optimo: [" << m_problem->getOptime() << "]" << endl;

    sol = p1;
    fitness = eval1;
}