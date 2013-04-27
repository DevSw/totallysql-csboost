using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class hypergeometric_distribution : distribution
    {
        ulong m_n, m_N, m_r;

        public hypergeometric_distribution(ulong r, ulong n, ulong N)
        {
            m_n = n;
            m_N = N;
            m_r = r;
            check_parameters();
        }
        public override bool discrete() { return true; }

        public ulong total() { return m_N; }

        public ulong defective() { return m_n; }

        public ulong sample_count() { return m_r; }

        public override bool symmetric()
        {
            return (m_n * 2 == m_N || m_r * 2 == m_N);
        }

        public override bool tail_right() { return false; }

        public override bool tail_left() { return false; }

        public override void setup_bars()
        {
            return;
        }

        public override void check_parameters()
        {
            if (m_r > m_N) throw new ArgumentException(string.Format("Hypergeometric distribution: parameter r must be <= N (got r={0:G}, N={1:G}).", m_r, m_N));
            if (m_n > m_N) throw new ArgumentException(string.Format("Hypergeometric distribution: parameter n must be <= N (got n={0:G}, N={1:G}).", m_n, m_N));
        }

        private void check_x(ulong x)
        {
            if(x < Math.Max(0, m_n + m_r - m_N)) throw new ArgumentException(string.Format("Random variable out of range: must be > 0 and > m + r - N but got {0:G}.", x));
            if(x > Math.Min(m_r, m_n)) throw new ArgumentException(string.Format("Random variable out of range: must be less than both n and r but got {0:G}.", x));
        }

        public override XMath.pair<double> range()
        {
            ulong l = m_n + m_r > m_N ? m_n + m_r - m_N : 0;
            ulong u = (ulong)Math.Min(m_r, m_n);
            return new XMath.pair<double>(l, u);
        }

        public override XMath.pair<double> support()
        {
            return range();
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            double result;
            if (m_N < XMath.max_factorial) result = hypergeometric_pdf_factorial_imp((ulong)x, m_r, m_n, m_N);
            else if (m_N <= XMath.max_prime_index) result = hypergeometric_pdf_prime_imp((ulong)x, m_r, m_n, m_N);
            else result = hypergeometric_pdf_lanczos_imp((ulong)x, m_r, m_n, m_N);

            if (result > 1) result = 1;
            if (result < 0) result = 0;

            return result;

        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (RHS) return find_pdf_inv(p, mode(), support().v2, false);
            else return find_pdf_inv(p, support().v1, mode(), true);
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            double result = 0;
            result = hypergeometric_cdf_imp((ulong)x, false);
            if (result > 1) result = 1;
            if (result < 0) result = 0;
            return result;
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            double result = 0;
            result = hypergeometric_cdf_imp((ulong)x, true);
            if (result > 1) result = 1;
            if (result < 0) result = 0;
            return result;
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            return hypergeometric_quantile_imp(p, 1 - p);
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            return hypergeometric_quantile_imp(1-q, q);
        }

        public override double mean()
        {
            return (double)m_n * (double)m_r / (double)m_N;
        }

        public override double variance()
        {
            return (double)m_r * ((double)m_n / (double)m_N) * (1.0 - (double)m_n / (double)m_N) * ((double)m_N - (double)m_r) / ((double)m_N - 1.0);
        }

        public override double mode()
        {
            return Math.Floor(((double)m_r + 1.0) * ((double)m_n + 1.0) / ((double)m_N + 2.0));
        }

        //Median supplied by base class

        public override double skewness()
        {
            double r = (double)m_r;
            double n = (double)m_n;
            double N = (double)m_N;
            return (N - 2 * n) * Math.Sqrt(N - 1) * (N - 2 * r) / (Math.Sqrt(n * r * (N - n) * (N - r)) * (N - 2));
        }

        public override double kurtosis_excess()
        {
            double r = (double)m_r;
            double n = (double)m_n;
            double N = (double)m_N;
            double t1 = N * N * (N - 1) / (r * (N - 2) * (N - 3) * (N - r));
            double t2 = (N * (N + 1) - 6 * N * (N - r)) / (n * (N - n))
               + 3 * r * (N - r) * (N + 6) / (N * N) - 6;
            return t1 * t2;
        }

        //Kurtosis supplied by base class

        public override double random()
        {
            return Math.Round(base.random());
        }

        private ulong hypergeometric_quantile_imp(double p, double q)
        {
            double result;
            double fudge_factor = 1 + XMath.epsilon * ((m_N <= XMath.max_prime_index) ? 50 : 2 * m_N);
            ulong bayse = (ulong)(Math.Max(0, (int)(m_n + m_r) - (int)(m_N)));
            ulong lim = (ulong)Math.Min(m_r, m_n);

            if (p <= 0.5)
            {
                ulong x = bayse;
                result = pdf(x);
                double diff = result;
                while (result < p)
                {
                    diff = (diff > XMath.min_value * 8) ? (double)(m_n - x) * (double)(m_r - x) * diff / ((x + 1.0) * (m_N + x + 1.0 - m_n - m_r)) : pdf(x);
                    if (result + diff / 2 > p)
                        break;
                    ++x;
                    result += diff;
                }
                return x;
            }
            else
            {
                ulong x = lim;
                result = 0;
                double diff = pdf(x);
                while (result + diff / 2 < q)
                {
                    result += diff;
                    diff = (diff > XMath.min_value * 8) ? x * (double)(m_N + x - m_n - m_r) * diff / ((1.0 + m_n - x) * (1.0 + m_r - x)) : pdf(x - 1.0);
                    --x;
                }
                return x;
            }
        }

        private double hypergeometric_cdf_imp(ulong x, bool invert)
        {
            double result = 0;
            double mode = Math.Floor((m_r + 1.0) * (m_n + 1.0) / (m_N + 2.0));
            if (x < mode)
            {
                result = pdf(x);
                double diff = result;
                ulong lower_limit = (ulong)(Math.Max(0, (int)(m_n + m_r) - (int)(m_N)));
                while (diff > (invert ? 1.0 : result) * XMath.epsilon)
                {
                    diff = (double)x * (double)((m_N + x) - m_n - m_r) * diff / ((1.0 + m_n - x) * (1.0 + m_r - x));
                    result += diff;
                    if (x == lower_limit)
                        break;
                    --x;
                }
            }
            else
            {
                invert = !invert;
                ulong upper_limit = (ulong)Math.Min(m_r, m_n);
                if (x != upper_limit)
                {
                    ++x;
                    result = pdf(x);
                    double diff = result;
                    while ((x <= upper_limit) && (diff > (invert ? 1.0 : result) * XMath.epsilon))
                    {
                        diff = (double)(m_n - x) * (double)(m_r - x) * diff / ((x + 1.0) * ((m_N + x + 1.0) - m_n - m_r));
                        result += diff;
                        ++x;
                    }
                }
            }
            if (invert)
                result = 1 - result;
            return result;
        }

        private double hypergeometric_pdf_factorial_imp(ulong x, ulong r, ulong n, ulong N)
        {
            double[] f = XMath.factorials();

            double result = f[n];
            double[] num = { f[r], f[N - n], f[N - r] };
            double[] denom = { f[N], f[x], f[n - x], f[r - x], f[N - n - r + x] };
            int i = 0;
            int j = 0;
            while ((i < 3) || (j < 5))
            {
                while ((j < 5) && ((result >= 1) || (i >= 3)))
                {
                    result /= denom[j];
                    ++j;
                }
                while ((i < 3) && ((result <= 1) || (j >= 5)))
                {
                    result *= num[i];
                    ++i;
                }
            }
            return result;
        }

        class hypergeometric_pdf_prime_loop_result_entry
        {
            public hypergeometric_pdf_prime_loop_result_entry(double v, hypergeometric_pdf_prime_loop_result_entry nxt)
            {
                value = v;
                next = nxt;
            }
            public double value;
            public hypergeometric_pdf_prime_loop_result_entry next;
        }

        class hypergeometric_pdf_prime_loop_data
        {
            public hypergeometric_pdf_prime_loop_data(ulong x_, ulong r_, ulong n_, ulong N_, int pi, ulong cp)
            {
                x = x_;
                r = r_;
                n = n_;
                N = N_;
                prime_index = pi;
                current_prime = cp;
            }
            public ulong x;
            public ulong r;
            public ulong n;
            public ulong N;
            public int prime_index;
            public ulong current_prime;
        }

        private double hypergeometric_pdf_prime_imp(ulong x, ulong r, ulong n, ulong N)
        {
            hypergeometric_pdf_prime_loop_result_entry result = new hypergeometric_pdf_prime_loop_result_entry(1, null);
            hypergeometric_pdf_prime_loop_data data = new hypergeometric_pdf_prime_loop_data(x, r, n, N, 0, XMath.prime(0));
            return hypergeometric_pdf_prime_loop_imp(data, result);
        }

        private double hypergeometric_pdf_prime_loop_imp(hypergeometric_pdf_prime_loop_data data, hypergeometric_pdf_prime_loop_result_entry result)
        {
            while (data.current_prime <= data.N)
            {
                ulong bayse = data.current_prime;
                int prime_powers = 0;
                while (bayse <= data.N)
                {
                    prime_powers += (int)(data.n / bayse);
                    prime_powers += (int)(data.r / bayse);
                    prime_powers += (int)((data.N - data.n) / bayse);
                    prime_powers += (int)((data.N - data.r) / bayse);
                    prime_powers -= (int)(data.N / bayse);
                    prime_powers -= (int)(data.x / bayse);
                    prime_powers -= (int)((data.n - data.x) / bayse);
                    prime_powers -= (int)((data.r - data.x) / bayse);
                    prime_powers -= (int)((data.N - data.n - data.r + data.x) / bayse);
                    bayse *= data.current_prime;
                }
                if (prime_powers != 0)
                {
                    double p = Math.Pow(data.current_prime, prime_powers);
                    if ((p > 1) && (double.MaxValue / p < result.value))
                    {
                        //
                        // The next calculation would overflow, use recursion
                        // to sidestep the issue:
                        //
                        hypergeometric_pdf_prime_loop_result_entry t = new hypergeometric_pdf_prime_loop_result_entry(p, result);
                        data.current_prime = XMath.prime(++data.prime_index);
                        return hypergeometric_pdf_prime_loop_imp(data, t);
                    }
                    if ((p < 1) && (ulong.MinValue / p > result.value))
                    {
                        //
                        // The next calculation would underflow, use recursion
                        // to sidestep the issue:
                        //
                        hypergeometric_pdf_prime_loop_result_entry t = new hypergeometric_pdf_prime_loop_result_entry(p, result);
                        data.current_prime = XMath.prime(++data.prime_index);
                        return hypergeometric_pdf_prime_loop_imp(data, t);
                    }
                    result.value *= p;
                }
                data.current_prime = XMath.prime(++data.prime_index);
            }
            //
            // When we get to here we have run out of prime factors,
            // the overall result is the product of all the partial
            // results we have accumulated on the stack so far, these
            // are in a linked list starting with "data.head" and ending
            // with "result".
            //
            // All that remains is to multiply them together, taking
            // care not to overflow or underflow.
            //
            // Enumerate partial results >= 1 in variable i
            // and partial results < 1 in variable j:
            //
            hypergeometric_pdf_prime_loop_result_entry i, j;
            i = result;
            while (i != null && i.value < 1) i = i.next;
            j = result;
            while (j != null && j.value >= 1) j = j.next;

            double prod = 1;

            while (i != null || j != null)
            {
                while (i != null && ((prod <= 1) || (j == null)))
                {
                    prod *= i.value;
                    i = i.next;
                    while (i != null && i.value < 1) i = i.next;
                }
                while (j != null && ((prod >= 1) || (i == null)))
                {
                    prod *= j.value;
                    j = j.next;
                    while (j != null && j.value >= 1) j = j.next;
                }
            }

            return prod;
        }

        struct sort_functor : IComparer<int>
        {
            List<double> m_exponents;

            public sort_functor(List<double> exponents)
            {
                m_exponents = exponents;
            }

            public int Compare(int i, int j)
            {
                return m_exponents[i].CompareTo(m_exponents[j]);
            }
        }

        static void bubble_down_one(List<int> list, int start, int end, sort_functor f)
        {
            int next = start;
            int first = start;
            ++next;
            while ((next != end) && (f.Compare(first, next) <= 0))
            {
                int hold = list[first];
                list[first] = list[next];
                list[next] = hold;
                ++first;
                ++next;
            }
        }

        double hypergeometric_pdf_lanczos_imp(ulong x, ulong r, ulong n, ulong N)
{

   List<double> bases = new List<double> (new double[]{
      n + XMath.lanczos_g + 0.5,
      r + XMath.lanczos_g + 0.5,
      N - n + XMath.lanczos_g + 0.5,
      N - r + XMath.lanczos_g + 0.5,
      1.0 / (N + XMath.lanczos_g + 0.5),
      1.0 / (x + XMath.lanczos_g + 0.5),
      1.0 / (n - x + XMath.lanczos_g + 0.5),
      1.0 / (r - x + XMath.lanczos_g + 0.5),
      1.0 / (N - n - r + x + XMath.lanczos_g + 0.5)
   });
   List<double> exponents = new List<double>( new double[] {
      n + 0.5,
      r + 0.5,
      N - n + 0.5,
      N - r + 0.5,
      N + 0.5,
      x + 0.5,
      n - x + 0.5,
      r - x + 0.5,
      N - n - r + x + 0.5
   });
   List<int> base_e_factors = new List<int>( new int[] {-1, -1, -1, -1, 1, 1, 1, 1, 1});
   List<int> sorted_indexes = new List<int>( new int[] {0, 1, 2, 3, 4, 5, 6, 7, 8});
            sort_functor f = new sort_functor(exponents);
            sorted_indexes.Sort(0, 9, f);

   do{
      exponents[sorted_indexes[0]] -= exponents[sorted_indexes[1]];
      bases[sorted_indexes[1]] *= bases[sorted_indexes[0]];
      if((bases[sorted_indexes[1]] < XMath.min_value) && (exponents[sorted_indexes[1]] != 0)) return 0;
      base_e_factors[sorted_indexes[1]] += base_e_factors[sorted_indexes[0]];
      bubble_down_one(sorted_indexes, 0, 9, f);
   }while(exponents[sorted_indexes[1]] > 1);

   //
   // Combine equal powers:
   //
   int j = 8;
   while(exponents[sorted_indexes[j]] == 0) --j;
   while(j != 0)
   {
      while(j != 0 && (exponents[sorted_indexes[j-1]] == exponents[sorted_indexes[j]]))
      {
         bases[sorted_indexes[j-1]] *= bases[sorted_indexes[j]];
         exponents[sorted_indexes[j]] = 0;
         base_e_factors[sorted_indexes[j-1]] += base_e_factors[sorted_indexes[j]];
         bubble_down_one(sorted_indexes, j, 9, f);
         --j;
      }
      --j;
   }

   double result;
   result = Math.Pow(bases[sorted_indexes[0]] * Math.Exp((double)(base_e_factors[sorted_indexes[0]])), exponents[sorted_indexes[0]]);
   for(int i = 1; (i < 9) && (exponents[sorted_indexes[i]] > 0); ++i)
   {
      if(result < XMath.min_value) return 0; // short circuit further evaluation
      if(exponents[sorted_indexes[i]] == 1)
         result *= bases[sorted_indexes[i]] * Math.Exp((double)(base_e_factors[sorted_indexes[i]]));
      else if(exponents[sorted_indexes[i]] == 0.5)
         result *= Math.Sqrt(bases[sorted_indexes[i]] * Math.Exp((double)(base_e_factors[sorted_indexes[i]])));
      else
         result *= Math.Pow(bases[sorted_indexes[i]] * Math.Exp((double)(base_e_factors[sorted_indexes[i]])), exponents[sorted_indexes[i]]);
   
   }

   result *= XMath.lanczos_sum_expG_scaled((double)(n + 1))
      * XMath.lanczos_sum_expG_scaled((double)(r + 1))
      * XMath.lanczos_sum_expG_scaled((double)(N - n + 1))
      * XMath.lanczos_sum_expG_scaled((double)(N - r + 1))
      / 
      ( XMath.lanczos_sum_expG_scaled((double)(N + 1))
         * XMath.lanczos_sum_expG_scaled((double)(x + 1))
         * XMath.lanczos_sum_expG_scaled((double)(n - x + 1))
         * XMath.lanczos_sum_expG_scaled((double)(r - x + 1))
         * XMath.lanczos_sum_expG_scaled((double)(N - n - r + x + 1)));
   
   return result;
}
    }
}
