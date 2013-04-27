using System;
using System.Collections.Generic;
using System.Text;

namespace CSBoost
{
    partial class XMath
    {

        public static double expint(uint n, double z)
        {
            return expint_imp(n, z);
        }

        public static double expint(double z)
        {
            return expint_i_imp(z);
        }

        #region implementation

        private static double expint_1_rational(double z)
        {
           double result;
           if(z <= 1)
           {
              // Maximum Deviation Found:                     2.006e-18
              // Expected Error Term:                         2.006e-18
              // Max error found at double precision:         2.760e-17
              double Y = 0.66373538970947265625;
              double[] P= {    
                 0.0865197248079397976498,
                 0.0320913665303559189999,
                 -0.245088216639761496153,
                 -0.0368031736257943745142,
                 -0.00399167106081113256961,
                 -0.000111507792921197858394
              };
              double[] Q= {    
                 1,
                 0.37091387659397013215,
                 0.056770677104207528384,
                 0.00427347600017103698101,
                 0.000131049900798434683324,
                 -0.528611029520217142048e-6
              };
              result = evaluate_polynomial(P, z) 
                 / evaluate_polynomial(Q, z);
              result += z - Math.Log(z) - Y;
           }
           else if(z < -double.Epsilon)
           {
              // Maximum Deviation Found (interpolated):      1.444e-17
              // Max error found at double precision:         3.119e-17
              double[] P= {    
                 -0.121013190657725568138e-18,
                 -0.999999999999998811143,
                 -43.3058660811817946037,
                 -724.581482791462469795,
                 -6046.8250112711035463,
                 -27182.6254466733970467,
                 -66598.2652345418633509,
                 -86273.1567711649528784,
                 -54844.4587226402067411,
                 -14751.4895786128450662,
                 -1185.45720315201027667
              };
              double[] Q= {    
                 1,
                 45.3058660811801465927,
                 809.193214954550328455,
                 7417.37624454689546708,
                 38129.5594484818471461,
                 113057.05869159631492,
                 192104.047790227984431,
                 180329.498380501819718,
                 86722.3403467334749201,
                 18455.4124737722049515,
                 1229.20784182403048905,
                 -0.776491285282330997549
              };
              double recip = 1 / z;
              result = 1 + evaluate_polynomial(P, recip)
                 / evaluate_polynomial(Q, recip);
              result *= Math.Exp(-z) * recip;
           }
           else
           {
              result = 0;
           }
           return result;
        }

        class expint_fraction : series<pair<double>>
        {
           public expint_fraction(uint n_, double z_)
           {
               b = n_ + z_;
               i = -1;
               n= n_;
           }

           public override pair<double> next()
           {
              pair<double> result = new pair<double>(-((i+1.0) * (n+i)), b);
              b += 2;
              ++i;
              return result;
           }
           double b;
           int i;
           uint n;
        }

        private static double expint_as_fraction(uint n, double z)
        {
           int max_iter = max_series_iterations;
           expint_fraction f = new expint_fraction(n, z);
           double result = continued_fraction_b(f, epsilon, ref max_iter);
           check_series_iterations((ulong)max_iter);
           result = Math.Exp(-z) / result;
           return result;
        }

        class expint_series : series<double>
        {
           public expint_series(uint k_, double z_, double x_k_, double denom_, double fact_) 
           {
               k = k_;
               z = z_;
               x_k = x_k_;
               denom = denom_;
               fact = fact_;
           }

            public override double next()
           {
              x_k *= -z;
              denom += 1;
              fact *= ++k;
              return x_k / (denom * fact);
           }

           uint k;
           double z;
           double x_k;
           double denom;
           double fact;
        }

        private static double expint_as_series(uint n, double z)
        {
           int max_iter = max_series_iterations;

           double result = 0;
           double x_k = -1;
           double denom = 1.0 - n;
           double fact = 1;
           uint k = 0;
           for(; k < n - 1;)
           {
              result += x_k / (denom * fact);
              denom += 1;
              x_k *= -z;
              fact *= ++k;
           }
           result += Math.Pow(-z, n - 1) * (digamma(n) - Math.Log(z)) / fact;

           expint_series s = new expint_series(k, z, x_k, denom, fact);
           result = sum_series(s, epsilon, ref max_iter, result);
           check_series_iterations((ulong)max_iter);
           return result;
        }

        private static double expint_imp(uint n, double z)
        {
           if(z < 0) throw new Exception(string.Format("Function requires z >= 0 but got {0:G}.", z));
           if(z == 0) throw new OverflowException();

           double result;

           bool f;
           if(n < 3)
           {
              f = z < 0.5;
           }
           else
           {
              f = z < ((n - 2.0) / (n - 1.0));
           }
           if(n == 0)
              result = Math.Exp(-z) / z;
           else if(n == 1) result = expint_1_rational(z);
           else if(f) result = expint_as_series(n, z);
           else result = expint_as_fraction(n, z);

           return result;
        }

        class expint_i_series : series<double>
        {
           public expint_i_series(double z_)
           {
               k = 0;
               z_k = 1;
               z = z_;
           }
           public override double next ()
           {
              z_k *= z / ++k;
              return z_k / k;
           }
           uint k;
           double z_k;
           double z;
        };

        public static double expint_i_as_series(double z)
        {
           double result = Math.Log(z); // (Math.Log(z) - Math.Log(1 / z)) / 2;
           result += euler;
           expint_i_series s = new expint_i_series(z);
           int max_iter = max_series_iterations;
           result = sum_series(s, epsilon, ref max_iter, result);
           check_series_iterations((ulong)max_iter);
           return result;
        }

        public static double expint_i_imp(double z)
        {
           if(z < 0)
              return -expint_imp(1, -z);
           if(z == 0) throw new OverflowException();

           double result;

           if(z <= 6)
           {
              // Maximum Deviation Found:                     2.852e-18
              // Expected Error Term:                         2.852e-18
              // Max Error found at double precision =        Poly: 2.636335e-16   Cheb: 4.187027e-16
              double[] P= {    
                 2.98677224343598593013,
                 0.356343618769377415068,
                 0.780836076283730801839,
                 0.114670926327032002811,
                 0.0499434773576515260534,
                 0.00726224593341228159561,
                 0.00115478237227804306827,
                 0.000116419523609765200999,
                 0.798296365679269702435e-5,
                 0.2777056254402008721e-6
              };
              double[] Q= {    
                 1,
                 -1.17090412365413911947,
                 0.62215109846016746276,
                 -0.195114782069495403315,
                 0.0391523431392967238166,
                 -0.00504800158663705747345,
                 0.000389034007436065401822,
                 -0.138972589601781706598e-4
              };

              double r1 = 1677624236387711.0 / 4503599627370496.0;
              double r2 = 0.131401834143860282009280387409357165515556574352422001206362e-16;
              double r = 0.372507410781366634461991866580119133535689497771654051555657435242200120636201854384926049951548942392;
              double t = (z / 3) - 1;
              result = evaluate_polynomial(P, t) 
                 / evaluate_polynomial(Q, t);
              t = (z - r1) - r2;
              result *= t;
              if(Math.Abs(t) < 0.1)
              {
                 result += log1p(t / r);
              }
              else
              {
                 result += Math.Log(z / r);
              }
           }
           else if (z <= 10)
           {
              // Maximum Deviation Found:                     6.546e-17
              // Expected Error Term:                         6.546e-17
              // Max Error found at double precision =        Poly: 6.890169e-17   Cheb: 6.772128e-17
              double Y = 1.158985137939453125;
              double[] P= {    
                 0.00139324086199402804173,
                 -0.0349921221823888744966,
                 -0.0264095520754134848538,
                 -0.00761224003005476438412,
                 -0.00247496209592143627977,
                 -0.000374885917942100256775,
                 -0.554086272024881826253e-4,
                 -0.396487648924804510056e-5
              };
              double[] Q= {    
                 1,
                 0.744625566823272107711,
                 0.329061095011767059236,
                 0.100128624977313872323,
                 0.0223851099128506347278,
                 0.00365334190742316650106,
                 0.000402453408512476836472,
                 0.263649630720255691787e-4
              };
              double t = z / 2 - 4;
              result = Y + evaluate_polynomial(P, t)
                 / evaluate_polynomial(Q, t);
              result *= Math.Exp(z) / z;
              result += z;
           }
           else if(z <= 20)
           {
              // Maximum Deviation Found:                     1.843e-17
              // Expected Error Term:                         -1.842e-17
              // Max Error found at double precision =        Poly: 4.375868e-17   Cheb: 5.860967e-17

              double Y = 1.0869731903076171875;
              double[] P= {    
                 -0.00893891094356945667451,
                 -0.0484607730127134045806,
                 -0.0652810444222236895772,
                 -0.0478447572647309671455,
                 -0.0226059218923777094596,
                 -0.00720603636917482065907,
                 -0.00155941947035972031334,
                 -0.000209750022660200888349,
                 -0.138652200349182596186e-4
              };
              double[] Q= {    
                 1,
                 1.97017214039061194971,
                 1.86232465043073157508,
                 1.09601437090337519977,
                 0.438873285773088870812,
                 0.122537731979686102756,
                 0.0233458478275769288159,
                 0.00278170769163303669021,
                 0.000159150281166108755531
              };
              double t = z / 5 - 3;
              result = Y + evaluate_polynomial(P, t)
                 / evaluate_polynomial(Q, t);
              result *= Math.Exp(z) / z;
              result += z;
           }
           else if(z <= 40)
           {
              // Maximum Deviation Found:                     5.102e-18
              // Expected Error Term:                         5.101e-18
              // Max Error found at double precision =        Poly: 1.441088e-16   Cheb: 1.864792e-16


              double Y = 1.03937530517578125;
              double[] P= {    
                 -0.00356165148914447597995,
                 -0.0229930320357982333406,
                 -0.0449814350482277917716,
                 -0.0453759383048193402336,
                 -0.0272050837209380717069,
                 -0.00994403059883350813295,
                 -0.00207592267812291726961,
                 -0.000192178045857733706044,
                 -0.113161784705911400295e-9
              };
              double[] Q= {    
                 1,
                 2.84354408840148561131,
                 3.6599610090072393012,
                 2.75088464344293083595,
                 1.2985244073998398643,
                 0.383213198510794507409,
                 0.0651165455496281337831,
                 0.00488071077519227853585
              };
              double t = z / 10 - 3;
              result = Y + evaluate_polynomial(P, t)
                 / evaluate_polynomial(Q, t);
              result *= Math.Exp(z) / z;
              result += z;
           }
           else
           {
              // Max Error found at double precision =        3.381886e-17
              double exp40 = 2.35385266837019985407899910749034804508871617254555467236651e17;
              double Y= 1.013065338134765625;
              double[] P= {    
                 -0.0130653381347656243849,
                 0.19029710559486576682,
                 94.7365094537197236011,
                 -2516.35323679844256203,
                 18932.0850014925993025,
                 -38703.1431362056714134
              };
              double[] Q= {    
                 1,
                 61.9733592849439884145,
                 -2354.56211323420194283,
                 22329.1459489893079041,
                 -70126.245140396567133,
                 54738.2833147775537106,
                 8297.16296356518409347
              };
              double t = 1 / z;
              result = Y + evaluate_polynomial(P, t)
                 / evaluate_polynomial(Q, t);
              if(z < 41)
                 result *= Math.Exp(z) / z;
              else
              {
                 // Avoid premature overflow if we can:
                 t = z - 40;
                 if(t > log_max_value) throw new OverflowException();
                 else
                 {
                    result *= Math.Exp(z - 40) / z;
                    if(result > double.MaxValue / exp40) throw new OverflowException();
                    else result *= exp40;
                 }
              }
              result += z;
           }
           return result;
        }

        #endregion
    }
}
