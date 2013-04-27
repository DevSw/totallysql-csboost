using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost
{
    partial class XMath
    {
        public static double zeta(double s, double sc)
        {
            if (s == 1) throw new ArgumentException("Zeta function is discontinuous at s = 1");
            double result;
            if (s == 0) result = -0.5;
            else if (s < 0)
            {
                swap(ref s, ref sc);
                if (Math.Floor(sc / 2) == sc / 2) result = 0;
                else
                {
                    result = sin_pi(0.5 * sc) * 2 * Math.Pow(2 * Math.PI, -s) * gamma(s) * zeta(s, sc);
                }
            }
            else
            {
                result = zeta_imp_prec(s, sc);
            }
            return result;
        }

        private static double zeta_imp_prec(double s, double sc)
        {
            double result;
            if (s < 1)
            {
                // Rational Approximation
                // Maximum Deviation Found:                     2.020e-18
                // Expected Error Term:                         -2.020e-18
                // Max error found at double precision:         3.994987e-17
                double[] P = 
                {    
                    0.24339294433593750202,
                    -0.49092470516353571651,
                    0.0557616214776046784287,
                    -0.00320912498879085894856,
                    0.000451534528645796438704,
                    -0.933241270357061460782e-5
                };

                double[] Q = 
                {
                    1,
                    -0.279960334310344432495,
                    0.0419676223309986037706,
                    -0.00413421406552171059003,
                    0.00024978985622317935355,
                    -0.101855788418564031874e-4
                };

                result = evaluate_polynomial(P, sc) / evaluate_polynomial(Q, sc);
                result -= 1.2433929443359375;
                result += (sc);
                result /= (sc);
            }
            else if (s <= 2)
            {
                // Maximum Deviation Found:        9.007e-20
                // Expected Error Term:            9.007e-20
                double[] P = 
                {
                    0.577215664901532860516,
                    0.243210646940107164097,
                    0.0417364673988216497593,
                    0.00390252087072843288378,
                    0.000249606367151877175456,
                    0.110108440976732897969e-4
                };

                double[] Q = 
                {
                    1,
                    0.295201277126631761737,
                    0.043460910607305495864,
                    0.00434930582085826330659,
                    0.000255784226140488490982,
                    0.10991819782396112081e-4,
                };

                result = evaluate_polynomial(P, -sc) / evaluate_polynomial(Q, -sc);
                result += 1 / (-sc);
            }
            else if (s <= 4)
            {
                // Maximum Deviation Found:          5.946e-22
                // Expected Error Term:              -5.946e-22
                double Y = 0.6986598968505859375;
                double[] P = 
                {
                    -0.0537258300023595030676,
                    0.0445163473292365591906,
                    0.0128677673534519952905,
                    0.00097541770457391752726,
                    0.769875101573654070925e-4,
                    0.328032510000383084155e-5,
                };
                double[] Q = 
                {
                    1,
                    0.33383194553034051422,
                    0.0487798431291407621462,
                    0.00479039708573558490716,
                    0.000270776703956336357707,
                    0.106951867532057341359e-4,
                    0.236276623974978646399e-7,
                };
                result = evaluate_polynomial(P, s - 2) / evaluate_polynomial(Q, s - 2);
                result += Y + 1 / (-sc);
            }
            else if (s <= 7)
            {
                // Maximum Deviation Found:                     2.955e-17
                // Expected Error Term:                         2.955e-17
                // Max error found at double precision:         2.009135e-16

                double[] P = 
                {
                    -2.49710190602259410021,
                    -2.60013301809475665334,
                    -0.939260435377109939261,
                    -0.138448617995741530935,
                    -0.00701721240549802377623,
                    -0.229257310594893932383e-4,
                };
                double[] Q = 
                {    
                    1,
                    0.706039025937745133628,
                    0.15739599649558626358,
                    0.0106117950976845084417,
                    -0.36910273311764618902e-4,
                    0.493409563927590008943e-5,
                    -0.234055487025287216506e-6,
                    0.718833729365459760664e-8,
                    -0.1129200113474947419e-9,
                };
                result = evaluate_polynomial(P, s - 4) / evaluate_polynomial(Q, s - 4);
                result = 1 + Math.Exp(result);
            }
            else if (s < 15)
            {
                // Maximum Deviation Found:                     7.117e-16
                // Expected Error Term:                         7.117e-16
                // Max error found at double precision:         9.387771e-16
                double[] P = 
                {    
                    -4.78558028495135619286,
                    -1.89197364881972536382,
                    -0.211407134874412820099,
                    -0.000189204758260076688518,
                    0.00115140923889178742086,
                    0.639949204213164496988e-4,
                    0.139348932445324888343e-5,
                };
                double[] Q = 
                {    
                    1,
                    0.244345337378188557777,
                    0.00873370754492288653669,
                    -0.00117592765334434471562,
                    -0.743743682899933180415e-4,
                    -0.21750464515767984778e-5,
                    0.471001264003076486547e-8,
                    -0.833378440625385520576e-10,
                    0.699841545204845636531e-12,
                };
                result = evaluate_polynomial(P, s - 7) / evaluate_polynomial(Q, s - 7);
                result = 1 + Math.Exp(result);
            }
            else if (s < 36)
            {
                // Max error in interpolated form:             1.668e-17
                // Max error found at long double precision:   1.669714e-17
                double[] P = 
                {    
                    -10.3948950573308896825,
                    -2.85827219671106697179,
                    -0.347728266539245787271,
                    -0.0251156064655346341766,
                    -0.00119459173416968685689,
                    -0.382529323507967522614e-4,
                    -0.785523633796723466968e-6,
                    -0.821465709095465524192e-8,
                };
                double[] Q = 
                {    
                    1,
                    0.208196333572671890965,
                    0.0195687657317205033485,
                    0.00111079638102485921877,
                    0.408507746266039256231e-4,
                    0.955561123065693483991e-6,
                    0.118507153474022900583e-7,
                    0.222609483627352615142e-14,
                };
                result = evaluate_polynomial(P, s - 15) / evaluate_polynomial(Q, s - 15);
                result = 1 + Math.Exp(result);
            }
            else if (s < 56)
            {
                result = 1 + Math.Pow(2.0, -s);
            }
            else
            {
                result = 1;
            }
            return result;
        }

    }
}
