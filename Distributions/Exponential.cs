using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class exponential_distribution: distribution
    {
        double m_lambda;

        public exponential_distribution(double l)
        {
            m_lambda = l;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_lambda <= 0 || double.IsInfinity(m_lambda)) throw new ArgumentException(string.Format("Lamda argument must be a finite number > 0 (got {0:G}).", m_lambda));
        }

        public override bool discrete() { return false; }

        public double lambda() { return m_lambda; }

        public override bool LHS() { return false; }

        public override bool tail_left() { return false; }

        public override bool strictly_decreasing()
        {
            return true;
        }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            return m_lambda * Math.Exp(-m_lambda * x);
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            if (!RHS) throw new Exception("The Exponential Distribution has no LHS (so no LHS inverse PDF is available");
            base.pdf_inv(p, RHS);
            if (p == 0) return double.MaxValue;
            return -Math.Log(p / m_lambda) / m_lambda;
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            return -XMath.expm1(-x * m_lambda);
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            return Math.Exp(-x * m_lambda);
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            if (p == 0) return 0;
            if (p == 1) throw new OverflowException();
            return -XMath.log1p(-p) / m_lambda;
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            if (q == 1) return 0;
            if (q == 0) throw new OverflowException();
            return -Math.Log(q) / m_lambda;
        }

        public override double mean()
        {
            return 1.0 / m_lambda;
        }

        //variance supplied by base class

        public override double standard_deviation()
        {
            return 1.0 / m_lambda;
        }

        public override double mode()
        {
            return 0.0;
        }

        public override double median()
        {
            return Math.Log(2.0) / m_lambda;
        }

        public override double skewness()
        {
            return 2.0;
        }

        //kurtosis supplied by base clase

        public override double kurtosis_excess()
        {
            return 6.0;
        }
    }
}
