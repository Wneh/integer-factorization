#include <iostream>
#include <gmpxx.h>
#include <vector>
#include <time.h>

//g++ --std=c++0x main.cpp -lgmpxx -lgmp -o factoring.out

using namespace std;

vector<mpz_class> smallPrimes;
vector<mpz_class> factors;
int primetest = 25;
unsigned int numBrents = 100000;
// long totalBrents = 500000;
gmp_randstate_t r_state;
int timeToRun = 14;
time_t startTime;
time_t endTime;
int pollardLimit = 50;
const mpz_class ONE(1);


mpz_class fnuc(const mpz_class& y,const mpz_class& c,const mpz_class& mod){
    return (y*y+c)%mod;
}

void gcd(mpz_class& x, mpz_class& y, mpz_class& z){
    mpz_gcd(x.get_mpz_t(),y.get_mpz_t(),z.get_mpz_t());
}

bool isPrime(mpz_class& x){
    return mpz_probab_prime_p(x.get_mpz_t(),primetest);
}

void random(mpz_class& x, mpz_class& n){
    //mpz_class temp;
    //mpz_sqrt(temp.get_mpz_t(),n.get_mpz_t());
    mpz_urandomm(x.get_mpz_t(),r_state,n.get_mpz_t());
}

/*mpz_class& min(mpz_class& x, mpz_class& y){
    return (x < y ? x : y);
}*/

void abs1(mpz_class& ret, mpz_class& x){
    mpz_abs(ret.get_mpz_t(),x.get_mpz_t());
}

mpz_class pollard_rho(mpz_class& n){

    //for(int i = 0 ; i < smallPrimes.size() ; ++i){
        mpz_class x = 2;
        mpz_class y = 2;
        mpz_class c = 1;
        mpz_class d = 1;
        mpz_class temp;
        int i = 0;
        while(d == 1){
            //check so there still is time (will run all time)
            x = fnuc(x,c,n);
            y = fnuc(fnuc(y,c,n), c, n);
            temp = x-y;
            gcd(d,n,temp); //d is return
            i++;
            if(i > pollardLimit){
                return 1;
            }
            
        }
        if(d != n){
            //cout << d << endl;
            return d;
        }
    //}
    return 1;
}

mpz_class bigPrimeFactorer(mpz_class& n){
    mpz_class tmp(211);
    for(int j = 0 ; j < 500 ; ++j){
        if( n % tmp == 0 ){
            return tmp;
        }
           mpz_nextprime(tmp.get_mpz_t(),tmp.get_mpz_t());
    }
    return 1;   
}

void ApowK(mpz_class& n, mpz_class& x){
    if(isPrime(x)){
        //cout << "APOWK!" << endl;
        while( n % x == 0 ){
            //cout << "Another ones goes..." << endl;
            factors.push_back(x);
            n = n / x;
        }
    }
    //cout << "left: " << n << endl;
}

bool testbit(mpz_class a, unsigned long int i) {
    return (i == mpz_scan1(a.get_mpz_t(), i))? 1 : 0; 
}


mpz_class william(mpz_class& n){

    mpz_class d, count, vs, vb,tempvar,p;
    int cc, i, k;
    double bp;

    unsigned int m = mpz_sizeinbase(n.get_mpz_t(),2); //random

    if(m > 88){
        return 1;
    }


    cc = 10;
    count = 0;
    p = 4;

    for(int index = 0; index < numBrents; ++index){
        if(time(NULL) > endTime){
        //if(totalBrents < 0){
            return 1;
        }
        //cout << "Outer" << endl;
        tempvar = p - 2;
        gcd(d,tempvar,n);
        if(1 < d && d < n){
            if(isPrime(d)){
                return d;
            }
            else{
                william(d);
            }
        }
        for(i = 1; i <= cc; i++){
            //cout << "Inner" << endl;
            count++;
            vs = p;
            vb = (p*p - 2)%n;
            //bp =  log2(count);
            //bp = logn(count, 2);
            bp = mpz_sizeinbase(count.get_mpz_t(),2);
            for(k = bp; k > -1; k--){
                if(testbit(count,k)){
                    vs = (vs*vb - p)%n;
                    vb = (vb * vb -2)%n;
                }
                else{
                    vb = (vs * vb - p)%n;
                    vs = (vs *vs -2)%n;
                }
            }
            p = vs;
        }
    }

    return 1;

}



mpz_class pollard_brent(mpz_class& n){
        mpz_class y; //random

        random(y,n);
        //mpz_class c(1);//c(y); //random
        //random(c,n);
        unsigned int m = mpz_sizeinbase(n.get_mpz_t(),2); //random

       if(m > 99){
        return 1;
       }
        //random(m,n);
        mpz_class g(1);
        unsigned int r = 1;
        mpz_class q(1);

        mpz_class x(y);
        mpz_class ys(y);

        mpz_class temp;
        mpz_class tempvar;
         while(g == 1){
                x = y;
                for(int i = 0 ; i < r ; ++i){
                    y = ((y*y)+ONE)%n; //((y*y)%n+c)%n;
                }
                unsigned int k = 0;
                while(k < r && g == 1){
                    //ys = y;
                //tempvar = r-k;
                //abs1(tempvar,tempvar);
                    for(int i = 0 ; i < min(m,r-k); ++i){
                        y = ((y*y)+ONE)%n;
                        q = (q*(x-y))%n;
                    }
                    gcd(g, q, n);
                    k = k + m;
                    if(k > numBrents){
                        //total Brents -= k;
                        return 1;
                    }
                }
                r = r << 1;
            }
            /*if(g == n){
             int i = 0;
                while(true){
                    ys = ((ys*ys)+ONE)%n;
                    tempvar = x - ys;
                    //abs1(temp, tempvar);
                    gcd(g, tempvar ,n);
                i++;
                    if(g > 1 || i >= pollardLimit){
                        break;
                    }
                }
            }*/

        if(g != n){
            return g;
        }
        
       //return bigPrimeFactorer(g);
       return 1;

}



bool smallPrimeFactorer(mpz_class& n){
    for(int j = 0 ; j < smallPrimes.size() ; ++j){
        if( n % smallPrimes[j] == 0 ){
            factors.push_back(smallPrimes[j]);
            n = n/smallPrimes[j];
            --j;
        }
    }
    if(n == 1){
        return false;
    } else {
        return true;
    }
}

void factor(mpz_class& n){
    if(isPrime(n)){
        factors.push_back(n);
    }

    else {
        mpz_class n1;
        mpz_class n2;

        //n1 = pollard_brent(n);
        //n1 = pollard_rho(n);
        n1 = william(n);
        if(time(NULL) > endTime){
        //if(totalBrents < 0){
            n1 = 1;
        }
        if(n1 == 1){
            cout << "fail" << endl;
            factors.clear();
        } else {
            n2 = n / n1;
            factor(n1);
            ApowK(n2,n1);
        if(factors.size() != 0 && n2 != 1){
            factor(n2);
        }
        }

    }
}


int main()
{
    unsigned long int seed;

    seed = time(NULL);//time perhaps?

    gmp_randinit_default(r_state);
    gmp_randseed_ui(r_state, seed);

    smallPrimes.reserve(30);
    smallPrimes.push_back(2);
    smallPrimes.push_back(3);
    smallPrimes.push_back(5);
    smallPrimes.push_back(7);
    smallPrimes.push_back(11);
    smallPrimes.push_back(13);
    smallPrimes.push_back(17);
    smallPrimes.push_back(19);
    smallPrimes.push_back(23);
    smallPrimes.push_back(29); //10
    smallPrimes.push_back(31);
    smallPrimes.push_back(37);
    smallPrimes.push_back(41);
    smallPrimes.push_back(43);
    smallPrimes.push_back(47);
    smallPrimes.push_back(53);
    smallPrimes.push_back(59);
    smallPrimes.push_back(61);
    smallPrimes.push_back(67);
    smallPrimes.push_back(71); //20
    smallPrimes.push_back(83);
    smallPrimes.push_back(89);
    smallPrimes.push_back(97);
    smallPrimes.push_back(101);
    smallPrimes.push_back(103);
    smallPrimes.push_back(107);
    smallPrimes.push_back(109);
    smallPrimes.push_back(113);
    smallPrimes.push_back(127);
    smallPrimes.push_back(131); //30
    /*smallPrimes.push_back(137);
    smallPrimes.push_back(139);
    smallPrimes.push_back(149);
    smallPrimes.push_back(151);
    smallPrimes.push_back(157);
    smallPrimes.push_back(163);
    smallPrimes.push_back(167);
    smallPrimes.push_back(173);
    smallPrimes.push_back(179);
    smallPrimes.push_back(181); //40
    smallPrimes.push_back(191);
    smallPrimes.push_back(193);
    smallPrimes.push_back(197);
    smallPrimes.push_back(199);*/


    factors.reserve(30);
    mpz_class test;
    startTime = time(NULL);
    endTime = startTime + timeToRun;

    for(int i = 0 ; i < 100 ; ++i){
        if(time(NULL) < endTime){
            //cout << totalBrents << " " << 0 << endl;
        //if(totalBrents > 0){
            cin >> test;
            if(smallPrimeFactorer(test)){
                factor(test);
            }
            for(int i = 0; i < factors.size() ; ++i){
                cout << factors[i] << endl;
            }
        } else {
            cout << "fail" <<endl;
        }
    factors.clear();
        cout << endl;
    }

    return 0;
}



