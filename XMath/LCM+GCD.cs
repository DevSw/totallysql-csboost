using System;
using System.Collections.Generic;
using System.Text;

namespace CSBoost
{
    public partial class XMath
    {
        public static long gcd(long a, long b)
        {
            while (true)
            {
                if (a == 0) return Math.Abs(b);
                b %= a;

                if (b == 0) return Math.Abs(a);
                a %= b;
            }
        }

        public static long lcm(long a, long b)
        {
            long temp = gcd( a, b );
            return ( temp != 0 ) ? ( a / temp * b ) : 0;
        }
    }
}
