using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class poisson_distribution: distribution
    {
        double m_l;

        public poisson_distribution(double mean)
        {
            m_l = mean;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_l <= 0 || double.IsInfinity(m_l)) throw new ArgumentException(string.Format("Mean must be a finite number > 0 (got {0:G}).", m_l));
        }

        public override bool discrete() { return true; }

        public override bool tail_right() { return false; }

        public override bool tail_left() { return false; }

        public override double mean()
        {
            return m_l;
        }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return range();
        }

        public override double mode()
        {
            return Math.Floor(mean());
        }

        // Median implemented via base class

        public override double variance()
        {
            return mean();
        }

        public override double skewness()
        {
            return 1.0 / Math.Sqrt(mean());
        }

        public override double kurtosis_excess()
        {
            return 1.0 / mean();
        }

        // kurtosis implemented via base class;

        public override double pdf(double x)
        {
            base.pdf(x);
            if (mean() == 0) return 0;
            if (x == 0) return Math.Exp(-mean());
            return XMath.gamma_p_derivative(x + 1, mean());
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (p == 0) return RHS ? double.MaxValue : 0;
            if (RHS) return find_pdf_inv(p, (int)mode(), (int)quantilec(double.Epsilon), false);
            return find_pdf_inv(p, 0, (int)mode(), true);
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            if (mean() == 0) return 0;
            if (x == 0) return Math.Exp(-mean());
            return XMath.gamma_q(x + 1, mean());
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            if (mean() == 0) return 1;
            if (x == 0) return -XMath.expm1(-mean());
            return XMath.gamma_p(x + 1, mean());
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            int max_iter = XMath.max_root_iterations;
            double guess, factor = 8;
            double z = mean();
            if (z < 1) guess = z;
            else guess = XMath.inverse_poisson_cornish_fisher(z, p, 1.0 - p);
            if (z > 5)
            {
                if (z > 1000) factor = 1.01;
                else if (z > 50) factor = 1.1;
                else if (guess > 10) factor = 1.25;
                else factor = 2;
                if (guess < 1.1) factor = 8;
            }
            return inverse_discrete_quantile(p, 1 - p, guess, factor, 1.0, max_iter);
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            int max_iter = XMath.max_root_iterations;
            double guess, factor = 8;
            double z = mean();
            if (z < 1) guess = z;
            else guess = XMath.inverse_poisson_cornish_fisher(z, 1.0 - q, q);
            if (z > 5)
            {
                if (z > 1000) factor = 1.01;
                else if (z > 50) factor = 1.1;
                else if (guess > 10) factor = 1.25;
                else factor = 2;
                if (guess < 1.1) factor = 8;
            }
            return inverse_discrete_quantile(1 - q, q, guess, factor, 1.0, max_iter);
        }

        public override double random()
        {
            return Math.Round(base.random());
        }
    }
}
