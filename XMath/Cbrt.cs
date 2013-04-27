using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost
{
    public static partial class XMath
    {
        public static double cbrt(double z)
        {
            double[] P = 
            { 
              0.37568269008611818,
              1.3304968705558024,
              -1.4897101632445036,
              1.2875573098219835,
              -0.6398703759826468,
              0.13584489959258635,
            };

            double[] correction = 
            {
              0.62996052494743658238360530363911,  // 2^-2/3
              0.79370052598409973737585281963615,  // 2^-1/3
              1,
              1.2599210498948731647672106072782,   // 2^1/3
              1.5874010519681994747517056392723,   // 2^2/3
            };

            if (double.IsInfinity(z)) throw new Exception(string.Format("Argument must be finite (got {0:G})", z));

            int i_exp, sign = 1;
            if (z < 0)
            {
                z = -z;
                sign = -sign;
            }
            if (z == 0)
                return 0;

            double guess = frexp(z, out i_exp);
            int original_i_exp = i_exp; // save for later
            guess = evaluate_polynomial(P, guess);
            int i_exp3 = i_exp / 3;

            if (Math.Abs(i_exp3) < 64)
            {
                if (i_exp3 > 0)
                    guess *= (ulong)1 << i_exp3;
                else
                    guess /= (ulong)1 << -i_exp3;
            }
            else
            {
                guess = ldexp(guess, i_exp3);
            }
            i_exp %= 3;
            guess *= correction[i_exp + 2];
            //
            // Now inline Halley iteration.
            // We do this here rather than calling tools::halley_iterate since we can
            // simplify the expressions algebraically, and don't need most of the error
            // checking of the boilerplate version as we know in advance that the function
            // is well behaved...
            //
            //
            // Epsilon calculation uses compile time arithmetic when it's available for type double,
            // otherwise uses ldexp to calculate at runtime:
            //
            double eps = ldexp(1.0, -2 - 53 / 3);
            double diff;

            if (original_i_exp < 1024 - 3)
            {
                //
                // Safe from overflow, use the fast method:
                //
                do
                {
                    double g3 = guess * guess * guess;
                    diff = (g3 + z + z) / (g3 + g3 + z);
                    guess *= diff;
                }
                while (Math.Abs(1 - diff) > eps);
            }
            else
            {
                //
                // Either we're ready to overflow, or we can't tell because numeric_limits isn't
                // available for type double:
                //
                do
                {
                    double g2 = guess * guess;
                    diff = (g2 - z / guess) / (2 * guess + z / g2);
                    guess -= diff;
                }
                while ((guess * eps) < Math.Abs(diff));
            }

            return sign * guess;
        }
    }
}
