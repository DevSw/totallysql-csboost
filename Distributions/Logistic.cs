using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class logistic_distribution : distribution
    {
        double m_location, m_scale;

        public logistic_distribution(double location, double scale)
        {
            m_location = location;
            m_scale = scale;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (double.IsInfinity(m_location)) throw new ArgumentException(string.Format("Location argument must be a finite number (got {0:G}).", m_location));
            if (m_scale <= 0 || double.IsInfinity(m_scale)) throw new ArgumentException(string.Format("Scale argument must be a finite number > 0 (got {0:G}).", m_scale));
        }

        public override bool discrete() { return false; }

        public override bool symmetric() { return true; }

        public double location() { return m_location; }

        public double scale() { return m_scale; }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(-double.MaxValue, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return range();
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            double exp_term = (m_location - x) / m_scale;
            if (Math.Abs(exp_term) > XMath.log_max_value) return 0;
            exp_term = Math.Exp(exp_term);
            if ((exp_term * m_scale > 1) && (exp_term > double.MaxValue / (m_scale * exp_term))) return 1 / (m_scale * exp_term);
            return (exp_term) / (m_scale * (1 + exp_term) * (1 + exp_term));
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (p == 0) return RHS ? double.MaxValue : -double.MaxValue;
            if (RHS) return find_pdf_inv(p, mode(), quantilec(XMath.epsilon), false);
            else return find_pdf_inv(p, quantile(XMath.epsilon), mode(), true);
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            if (double.IsNegativeInfinity(x)) return 0;
            if (double.IsPositiveInfinity(x)) return 1;
            double power = (m_location - x) / m_scale;
            if (power > XMath.log_max_value) return 0;
            if (power < -XMath.log_max_value) return 1;
            return 1 / (1 + Math.Exp(power));
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            if (double.IsNegativeInfinity(x)) return 1;
            if (double.IsPositiveInfinity(x)) return 0;
            double power = (x - m_location) / m_scale;
            if (power > XMath.log_max_value) return 0;
            if (power < -XMath.log_max_value) return 1;
            return 1 / (1 + Math.Exp(power));
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            if (p == 0 || p == 1) throw new OverflowException();
            return m_location - m_scale * Math.Log((1 - p) / p);
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            if (q == 0 || q == 1) throw new OverflowException();
            return m_location + m_scale * Math.Log((1 - q) / q);
        }

        public override double mean()
        {
            return m_location;
        }

        public override double variance()
        {
            return Math.PI * Math.PI * m_scale * m_scale / 3;
        }

        public override double mode()
        {
            return m_location;
        }

        public override double median()
        {
            return m_location;
        }

        public override double skewness()
        {
            return 0;
        }

        public override double kurtosis_excess()
        {
            return 6.0 / 5;
        }
    }
}
