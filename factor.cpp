#include <iostream>
#include <gmpxx.h>
#include <vector>
#include <time.h>

/* Compile with:
 * g++ --std=c++0x factor.cpp -lgmpxx -lgmp
 */


using namespace std;

vector<mpz_class> smallPrimes;
vector<mpz_class> factors;
vector<mpz_class> primes;
int primetest = 25;
unsigned int numBrents = 110000;
gmp_randstate_t r_state;
int timeToRun = 15;
time_t startTime;
time_t endTime;
int pollardLimit = 50;
const mpz_class ONE(1);

/* Gain next random for Pollard's rho. */
mpz_class fnuc(const mpz_class& y,const mpz_class& c,const mpz_class& mod){
    return (y*y+c)%mod;
}

/* Get GCD for y and z */
void gcd(mpz_class& x, mpz_class& y, mpz_class& z){
    mpz_gcd(x.get_mpz_t(),y.get_mpz_t(),z.get_mpz_t());
}

/* Milley-Rabin prime test */
bool isPrime(mpz_class& x){
    return mpz_probab_prime_p(x.get_mpz_t(),primetest);
}

/* Get random value between 0 and sqrt(n) */
void random(mpz_class& x, mpz_class& n){
    mpz_class temp;
    mpz_sqrt(temp.get_mpz_t(),n.get_mpz_t());
    mpz_urandomm(x.get_mpz_t(),r_state,temp.get_mpz_t());
}

/* Min between two mpz_class. */
mpz_class& min(mpz_class& x, mpz_class& y){
    return (x < y ? x : y);
}

/* Abs of two mpz_class. */
void abs1(mpz_class& ret, mpz_class& x){
    mpz_abs(ret.get_mpz_t(),x.get_mpz_t());
}

/* This pollard rho grants 68 points in kattis */
mpz_class pollard_rho(mpz_class& n){	
	mpz_class x = 2; //Should really be random but w/e
	mpz_class y = 2;
	mpz_class c = 1;//g++ --std=c++0x main.cpp -lgmpxx -lgmp -o factoring.out
	mpz_class d = 1;
	mpz_class temp;
	int i = 0;
	while(d == 1){
		x = fnuc(x,c,n);
		y = fnuc(fnuc(y,c,n), c, n);
		temp = x-y;
		gcd(d,n,temp);
		i++;
		if(i > pollardLimit){
			return 1;
		}
		
	}
	if(d != n){
		return d;
	}
	return 1;
}

/* Attempt factor with trial division of 500 primes. VERY SLOW */
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

/* Gains a list of every prime up to 200000. */
void sieveOfEratosthenes(){
        vector<bool> tmp;
        tmp.reserve(200000);
        
        /*for(unsigned int i = 0 ; i < 200000 ; ++i){ //This is to make sure the vector only contains true. However, this does not seem to be a problem.
                tmp.push_back(true);
        }*/
        for(unsigned int i = 2 ; i < 448 ; ++i){ // 448 ~ sqrt(200000)
                if(tmp[i]){
                        for(int j = i*i; j < 200000; j += i){ // i^2 + (iter*i)
                                tmp[j] = false;
                        }
                }
        }
        vector<mpz_class> ret;
        for(unsigned int i = 2 ; i < 200000 ; ++i){
                if(tmp[i]){
                        ret.push_back((mpz_class)i);
                }
        }
        primes = ret;
}

/* grant a log base q even though numbers are represented in base 2 */
void logq(mpz_class log, const mpz_class q,const mpz_class n){
        int logn = mpz_sizeinbase(n.get_mpz_t(),2);
        int logq = mpz_sizeinbase(q.get_mpz_t(),2);
        log = (logn/logq);
}

/* Pollard P-1 fatorization method. Used as a second (fast as limited to primes < 200000) method incase brent fails to gain a factor. */
mpz_class pollard_p1(mpz_class& n){
        mpz_class c(2);
        for(int i = 0 ; i < primes.size() ; ++i){
			mpz_class pp = primes[i];
			while(pp < 200000){
				mpz_powm(c.get_mpz_t(),c.get_mpz_t(),primes[i].get_mpz_t(),n.get_mpz_t());
				pp *= primes[i];
			}
		}
		mpz_class g;
		c--;
		gcd(g,c,n);
		if(g > 1 && g < n){
			return g;
		}
		return 1;
        
}

/* Pollard Brent factorization method. Pretty much the same as Pollard's rho though uses brent cycle detection. */
mpz_class pollard_brent(mpz_class& n){
        mpz_class y; //random
        random(y,n);

        unsigned int m = mpz_sizeinbase(n.get_mpz_t(),2);
        
        unsigned int brentLimit = numBrents;
        if(m > 99){
    	   brentLimit = brentLimit >> 3;
	    } else if(m > 98){
           brentLimit = brentLimit >> 2;
	    } else if(m > 95){
           brentLimit = brentLimit >> 1;   
	    }


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
                    y = ((y*y)+ONE)%n;
                }
                unsigned int k = 0;
                while(k < r && g == 1){
                    for(int i = 0 ; i < min(m,r-k); ++i){
                        y = ((y*y)+ONE)%n;
                        q = (q*(x-y))%n;
                    }
					gcd(g, q, n);
					k = k + m;
					if(k > brentLimit){
						return pollard_p1(n); //could call for pollard rho, pollard p-1, bigPrimeFactorer etc though this will require more time.
					}
				}
                r = r << 1;
            }
        if(g != n){
            return g;
        }
        
	   return pollard_p1(n); //could call for pollard rho, pollard p-1, bigPrimeFactorer etc though this will require more time.

}


/* trial division of small factors */
bool smallPrimeFactorer(mpz_class& n){
	for(int j = 0 ; j < smallPrimes.size() ; ++j){
		if( n % smallPrimes[j] == 0 ){
			factors.push_back(smallPrimes[j]);
			n = n/smallPrimes[j];
			--j;
		}
	}
	//Do we need to continue factor?
	if(n == 1){
		return false;
	} else {
		return true;
	}
}

void factor(mpz_class& n){
	
	//is n prime? then we're done.
    if(isPrime(n)){
        factors.push_back(n);
    }

    else {
        mpz_class n1;
        mpz_class n2;

        n1 = pollard_brent(n);
        //n1 = pollard_p12(n);
        
        //Still time to continue or should we abbandon results?
        if(time(NULL) > endTime){
            n1 = 1;
        }
        if(n1 == 1){ //unable to factor in time
            cout << "fail" << endl;
            factors.clear();
            
        } else { //continue factor problem
            n2 = n / n1;
            factor(n1);
		if(factors.size() != 0 && n2 != 1){ //in case of failure
		    factor(n2);
		}
        }

    }
}


int main()
{
    unsigned long int seed;

    seed = time(NULL);

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

	sieveOfEratosthenes();
    factors.reserve(30);
    mpz_class test;
    startTime = time(NULL);
    endTime = startTime + timeToRun;

    for(int i = 0 ; i < 100 ; ++i){
        if(time(NULL) < endTime){
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
