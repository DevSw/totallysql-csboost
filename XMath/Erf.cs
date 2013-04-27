using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost
{
    public static partial class XMath
    {
        public static double erf(double z)
        {
            return erf_imp(z, false);
        }

        public static double erfc(double z)
        {
            return erf_imp(z, true);
        }

        public static double erfc_inv(double z)
        {
            if ((z < 0) || (z > 2))
                throw new Exception(string.Format("Argument outside range [0,2] in inverse erfc function (got p={0:G}).", z));
            if (z == 0)
                throw new OverflowException();
            if (z == 2)
                throw new OverflowException();
            //
            // Normalise the input, so it's in the range [0,1], we will
            // negate the result if z is outside that range.  This is a simple
            // application of the erfc reflection formula: erfc(-z) = 2 - erfc(z)
            //
            double p, q, s;
            if (z > 1)
            {
                q = 2 - z;
                p = 1 - q;
                s = -1;
            }
            else
            {
                p = 1 - z;
                q = z;
                s = 1;
            }

            //
            // And get the result, negating where required:
            //
            return s * erf_inv_imp(p, q);
        }

        public static double erf_inv(double z)
        {
            if ((z < -1) || (z > 1))
                throw new Exception(string.Format("Argument outside range [-1, 1] in inverse erf function (got p={0:G}).", z));
            if (z == 1)
                throw new OverflowException();
            if (z == -1)
                throw new OverflowException();
            if (z == 0)
                return 0;
            //
            // Normalise the input, so it's in the range [0,1], we will
            // negate the result if z is outside that range.  This is a simple
            // application of the erf reflection formula: erf(-z) = -erf(z)
            //
            double p, q, s;
            if (z < 0)
            {
                p = -z;
                q = 1 - p;
                s = -1;
            }
            else
            {
                p = z;
                q = 1 - z;
                s = 1;
            }
            //
            // And get the result, negating where required:
            //
            return s * erf_inv_imp(p, q);
        }

        #region implementation

        private static double erf_imp(double z, bool invert)
        {
            if (z < 0)
            {
                if (!invert) return -erf_imp(-z, invert);
                else if (z < -0.5) return 2 - erf_imp(-z, invert);
                else return 1 + erf_imp(-z, false);
            }

            double result;

            //
            // Big bunch of selection statements now to pick
            // which implementation to use,
            // try to put most likely options first:
            //
            if (z < 0.5)
            {
                //
                // We're going to calculate erf:
                //
                if (z < 1e-10)
                {
                    if (z == 0)
                    {
                        result = 0.0;
                    }
                    else
                    {
                        result = z * 1.125 + z * 0.003379167095512573896158903121545171688;
                    }
                }
                else
                {
                    // Maximum Deviation Found:                     1.561e-17
                    // Expected Error Term:                         1.561e-17
                    // Maximum Relative Change in Control Points:   1.155e-04
                    // Max Error found at double precision =        2.961182e-17

                    const double Y = 1.044948577880859375;
                    double[] P = 
                 { 0.0834305892146531832907,        -0.338165134459360935041,       -0.0509990735146777432841,
                   -0.00772758345802133288487,      -0.000322780120964605683831
                 };

                    double[] Q = 
                 {    
                    1,
                    0.455004033050794024546,
                    0.0875222600142252549554,
                    0.00858571925074406212772,
                    0.000370900071787748000569
                 };

                    double zz = z * z;
                    result = z * (Y + evaluate_polynomial(P, zz) / evaluate_polynomial(Q, zz));
                }
            }
            else if (invert ? (z < 28) : (z < 5.8))
            {
                //
                // We'll be calculating erfc:
                //
                invert = !invert;
                if (z < 1.5)
                {
                    // Maximum Deviation Found:                     3.702e-17
                    // Expected Error Term:                         3.702e-17
                    // Maximum Relative Change in Control Points:   2.845e-04
                    // Max Error found at double precision =        4.841816e-17
                    const double Y = 0.405935764312744140625;
                    double[] P = 
                 {  -0.098090592216281240205,   0.178114665841120341155,    0.191003695796775433986,
                    0.0888900368967884466578,   0.0195049001251218801359,   0.00180424538297014223957,
                 };

                    double[] Q = 
                 {    
                    1,  1.84759070983002217845,    1.42628004845511324508,      0.578052804889902404909,    
                    0.12385097467900864233,        0.0113385233577001411017,    0.337511472483094676155e-5,
                 };

                    result = Y + evaluate_polynomial(P, z - 0.5) / evaluate_polynomial(Q, z - 0.5);
                    result *= Math.Exp(-z * z) / z;
                }
                else if (z < 2.5)
                {
                    // Max Error found at double precision =        6.599585e-18
                    // Maximum Deviation Found:                     3.909e-18
                    // Expected Error Term:                         3.909e-18
                    // Maximum Relative Change in Control Points:   9.886e-05
                    const double Y = 0.50672817230224609375;
                    double[] P = {    
                    -0.0243500476207698441272,
                    0.0386540375035707201728,
                    0.04394818964209516296,
                    0.0175679436311802092299,
                    0.00323962406290842133584,
                    0.000235839115596880717416,
                 };
                    double[] Q = {    
                    1,
                    1.53991494948552447182,
                    0.982403709157920235114,
                    0.325732924782444448493,
                    0.0563921837420478160373,
                    0.00410369723978904575884,
                 };
                    result = Y + evaluate_polynomial(P, z - 1.5) / evaluate_polynomial(Q, z - 1.5);
                    result *= Math.Exp(-z * z) / z;
                }
                else if (z < 4.5)
                {
                    // Maximum Deviation Found:                     1.512e-17
                    // Expected Error Term:                         1.512e-17
                    // Maximum Relative Change in Control Points:   2.222e-04
                    // Max Error found at double precision =        2.062515e-17
                    const double Y = 0.5405750274658203125;
                    double[] P = {    
                    0.00295276716530971662634,
                    0.0137384425896355332126,
                    0.00840807615555585383007,
                    0.00212825620914618649141,
                    0.000250269961544794627958,
                    0.113212406648847561139e-4,
                 };
                    double[] Q = {    
                    1,
                    1.04217814166938418171,
                    0.442597659481563127003,
                    0.0958492726301061423444,
                    0.0105982906484876531489,
                    0.000479411269521714493907,
                 };
                    result = Y + evaluate_polynomial(P, z - 3.5) / evaluate_polynomial(Q, z - 3.5);
                    result *= Math.Exp(-z * z) / z;
                }
                else
                {
                    // Max Error found at double precision =        2.997958e-17
                    // Maximum Deviation Found:                     2.860e-17
                    // Expected Error Term:                         2.859e-17
                    // Maximum Relative Change in Control Points:   1.357e-05
                    const double Y = 0.5579090118408203125;
                    double[] P = {    
                    0.00628057170626964891937,
                    0.0175389834052493308818,
                    -0.212652252872804219852,
                    -0.687717681153649930619,
                    -2.5518551727311523996,
                    -3.22729451764143718517,
                    -2.8175401114513378771,
                 };
                    double[] Q = {    
                    1,
                    2.79257750980575282228,
                    11.0567237927800161565,
                    15.930646027911794143,
                    22.9367376522880577224,
                    13.5064170191802889145,
                    5.48409182238641741584,
                 };
                    result = Y + evaluate_polynomial(P, 1 / z) / evaluate_polynomial(Q, 1 / z);
                    result *= Math.Exp(-z * z) / z;
                }
            }
            else
            {
                //
                // Any value of z larger than 28 will underflow to zero:
                //
                result = 0;
                invert = !invert;
            }

            if (invert)
            {
                result = 1 - result;
            }

            return result;
        }

        private static double erf_inv_imp(double p, double q)
        {

            double result = 0;

            if (p <= 0.5)
            {
                //
                // Evaluate inverse erf using the rational approximation:
                //
                // x = p(p+10)(Y+R(p))
                //
                // Where Y is a constant, and R(p) is optimised for a low
                // absolute error compared to |Y|.
                //
                // double: Max error found: 2.001849e-18
                // long double: Max error found: 1.017064e-20
                // Maximum Deviation Found (actual error term at infinite precision) 8.030e-21
                //
                double Y = 0.0891314744949340820313;
                double[] P = 
                {    
                    -0.000508781949658280665617,
                    -0.00836874819741736770379,
                    0.0334806625409744615033,
                    -0.0126926147662974029034,
                    -0.0365637971411762664006,
                    0.0219878681111168899165,
                    0.00822687874676915743155,
                    -0.00538772965071242932965
                };
                double[] Q = 
                {    
                    1.0,
                    -0.970005043303290640362,
                    -1.56574558234175846809,
                    1.56221558398423026363,
                    0.662328840472002992063,
                    -0.71228902341542847553,
                    -0.0527396382340099713954,
                    0.0795283687341571680018,
                    -0.00233393759374190016776,
                    0.000886216390456424707504
                };
                double g = p * (p + 10);
                double r = evaluate_polynomial(P, p) / evaluate_polynomial(Q, p);
                result = g * Y + g * r;
            }
            else if (q >= 0.25)
            {
                //
                // Rational approximation for 0.5 > q >= 0.25
                //
                // x = Math.Sqrt(-2*Math.Log(q)) / (Y + R(q))
                //
                // Where Y is a constant, and R(q) is optimised for a low
                // absolute error compared to Y.
                //
                // double : Max error found: 7.403372e-17
                // long double : Max error found: 6.084616e-20
                // Maximum Deviation Found (error term) 4.811e-20
                //
                double Y = 2.249481201171875;
                double[] P = 
                {    
                    -0.202433508355938759655,
                    0.105264680699391713268,
                    8.37050328343119927838,
                    17.6447298408374015486,
                    -18.8510648058714251895,
                    -44.6382324441786960818,
                    17.445385985570866523,
                    21.1294655448340526258,
                    -3.67192254707729348546
                };
                double[] Q = 
                {    
                    1.0,
                    6.24264124854247537712,
                    3.9713437953343869095,
                    -28.6608180499800029974,
                    -20.1432634680485188801,
                    48.5609213108739935468,
                    10.8268667355460159008,
                    -22.6436933413139721736,
                    1.72114765761200282724
                };
                double g = Math.Sqrt(-2 * Math.Log(q));
                double xs = q - 0.25;
                double r = evaluate_polynomial(P, xs) / evaluate_polynomial(Q, xs);
                result = g / (Y + r);
            }
            else
            {
                //
                // For q < 0.25 we have a series of rational approximations all
                // of the general form:
                //
                // let: x = Math.Sqrt(-Math.Log(q))
                //
                // Then the result is given by:
                //
                // x(Y+R(x-B))
                //
                // where Y is a constant, B is the lowest value of x for which 
                // the approximation is valid, and R(x-B) is optimised for a low
                // absolute error compared to Y.
                //
                // Note that almost all code will really go through the first
                // or maybe second approximation.  After than we're dealing with very
                // small input values indeed: 80 and 128 bit long double's go all the
                // way down to ~ 1e-5000 so the "tail" is rather long...
                //
                double x = Math.Sqrt(-Math.Log(q));
                if (x < 3)
                {
                    // Max error found: 1.089051e-20
                    double Y = 0.807220458984375;
                    double[] P = 
                    {    
                        -0.131102781679951906451,
                        -0.163794047193317060787,
                        0.117030156341995252019,
                        0.387079738972604337464,
                        0.337785538912035898924,
                        0.142869534408157156766,
                        0.0290157910005329060432,
                        0.00214558995388805277169,
                        -0.679465575181126350155e-6,
                        0.285225331782217055858e-7,
                        -0.681149956853776992068e-9
                    };
                    double[] Q = 
                    {    
                        1.0,
                        3.46625407242567245975,
                        5.38168345707006855425,
                        4.77846592945843778382,
                        2.59301921623620271374,
                        0.848854343457902036425,
                        0.152264338295331783612,
                        0.01105924229346489121
                    };
                    double xs = x - 1.125;
                    double R = evaluate_polynomial(P, xs) / evaluate_polynomial(Q, xs);
                    result = Y * x + R * x;
                }
                else if (x < 6)
                {
                    // Max error found: 8.389174e-21
                    double Y = 0.93995571136474609375;
                    double[] P = 
                    {    
                        -0.0350353787183177984712,
                        -0.00222426529213447927281,
                        0.0185573306514231072324,
                        0.00950804701325919603619,
                        0.00187123492819559223345,
                        0.000157544617424960554631,
                        0.460469890584317994083e-5,
                        -0.230404776911882601748e-9,
                        0.266339227425782031962e-11
                    };
                    double[] Q = 
                    {    
                        1.0,
                        1.3653349817554063097,
                        0.762059164553623404043,
                        0.220091105764131249824,
                        0.0341589143670947727934,
                        0.00263861676657015992959,
                        0.764675292302794483503e-4
                    };
                    double xs = x - 3;
                    double R = evaluate_polynomial(P, xs) / evaluate_polynomial(Q, xs);
                    result = Y * x + R * x;
                }
                else if (x < 18)
                {
                    // Max error found: 1.481312e-19
                    double Y = 0.98362827301025390625;
                    double[] P = 
                    {    
                        -0.0167431005076633737133,
                        -0.00112951438745580278863,
                        0.00105628862152492910091,
                        0.000209386317487588078668,
                        0.149624783758342370182e-4,
                        0.449696789927706453732e-6,
                        0.462596163522878599135e-8,
                        -0.281128735628831791805e-13,
                        0.99055709973310326855e-16
                    };
                    double[] Q = 
                    {    
                        1.0,
                        0.591429344886417493481,
                        0.138151865749083321638,
                        0.0160746087093676504695,
                        0.000964011807005165528527,
                        0.275335474764726041141e-4,
                        0.282243172016108031869e-6
                    };
                    double xs = x - 6;
                    double R = evaluate_polynomial(P, xs) / evaluate_polynomial(Q, xs);
                    result = Y * x + R * x;
                }
                else if (x < 44)
                {
                    // Max error found: 5.697761e-20
                    double Y = 0.99714565277099609375;
                    double[] P = 
                    {    
                        -0.0024978212791898131227,
                        -0.779190719229053954292e-5,
                        0.254723037413027451751e-4,
                        0.162397777342510920873e-5,
                        0.396341011304801168516e-7,
                        0.411632831190944208473e-9,
                        0.145596286718675035587e-11,
                        -0.116765012397184275695e-17
                    };
                    double[] Q = 
                    {    
                        1.0,
                        0.207123112214422517181,
                        0.0169410838120975906478,
                        0.000690538265622684595676,
                        0.145007359818232637924e-4,
                        0.144437756628144157666e-6,
                        0.509761276599778486139e-9
                    };
                    double xs = x - 18;
                    double R = evaluate_polynomial(P, xs) / evaluate_polynomial(Q, xs);
                    result = Y * x + R * x;
                }
                else
                {
                    // Max error found: 1.279746e-20
                    double Y = 0.99941349029541015625;
                    double[] P = 
                    {    
                        -0.000539042911019078575891,
                        -0.28398759004727721098e-6,
                        0.899465114892291446442e-6,
                        0.229345859265920864296e-7,
                        0.225561444863500149219e-9,
                        0.947846627503022684216e-12,
                        0.135880130108924861008e-14,
                        -0.348890393399948882918e-21
                    };
                    double[] Q = 
                    {    
                        1.0,
                        0.0845746234001899436914,
                        0.00282092984726264681981,
                        0.468292921940894236786e-4,
                        0.399968812193862100054e-6,
                        0.161809290887904476097e-8,
                        0.231558608310259605225e-11
                    };
                    double xs = x - 44;
                    double R = evaluate_polynomial(P, xs) / evaluate_polynomial(Q, xs);
                    result = Y * x + R * x;
                }
            }
            return result;
        }

        #endregion
    }
}
