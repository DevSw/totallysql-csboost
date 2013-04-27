using System;
using System.Collections.Generic;
using System.Text;

namespace CSBoost
{
    public struct Complex
    {
        private double r_, i_;

        public double Imaginary
        {
            get { return i_; }
            set { i_ = value; }
        }

        public double Real
        {
            get { return r_; }
            set { r_ = value; }
        }

        public double Magnitude
        {
            get { return Abs(this); }
            set { Real = value * Math.Cos(Phase); Imaginary = value * Math.Sin(Phase); }
        }

        public double Phase
        {
            get { return Arg(this); }
            set { Real = Magnitude * Math.Cos(value); Imaginary = Magnitude * Math.Sin(value); }
        }

        public Complex Conjugate()
        {
            return new Complex(Real, -Imaginary);
        }

        public Complex(double real, double imag)
        {
            r_ = real;
            i_ = imag;
        }

        public static Complex FromPolarCoordinates(double r, double theta)
        {
            return new Complex(r * Math.Cos(theta), r * Math.Sin(theta));
        }

        public static Complex operator +(Complex a, Complex b)
        {
            Complex c = new Complex(a.Real + b.Real, a.Imaginary + b.Imaginary);
            return c;
        }

        public static Complex operator +(double a, Complex b)
        {
            Complex c = new Complex(a + b.Real, b.Imaginary);
            return c;
        }

        public static Complex operator +(Complex a, double b)
        {
            Complex c = new Complex(a.Real + b, a.Imaginary);
            return c;
        }

        public static Complex operator -(Complex a, Complex b)
        {
            Complex c = new Complex(a.Real - b.Real, a.Imaginary - b.Imaginary);
            return c;
        }

        public static Complex operator -(Complex a, double b)
        {
            Complex c = new Complex(a.Real - b, a.Imaginary);
            return c;
        }

        public static Complex operator -(double a, Complex b)
        {
            Complex c = new Complex(a - b.Real, -b.Imaginary);
            return c;
        }

        public static Complex operator -(Complex c)
        {
            return 0 - c;
        }

        public static Complex operator *(Complex a, Complex b)
        {
            Complex c = new Complex(a.Real * b.Real - a.Imaginary * b.Imaginary, a.Real * b.Imaginary + a.Imaginary * b.Real);
            return c;
        }

        public static Complex operator *(Complex a, double b)
        {
            Complex c = new Complex(a.Real * b, a.Imaginary * b);
            return c;
        }

        public static Complex operator *(double a, Complex b)
        {
            Complex c = new Complex(a * b.Real, a * b.Imaginary);
            return c;
        }

        public static Complex operator /(Complex a, Complex b)
        {
            double denom = b.Real * b.Real + b.Imaginary * b.Imaginary;
            Complex c = new Complex((a.Real * b.Real + a.Imaginary * b.Imaginary) / denom, (a.Imaginary * b.Real + -a.Real * b.Imaginary) / denom);
            return c;
        }

        public static Complex operator /(Complex a, double b)
        {
            Complex c = new Complex(a.Real / b, a.Imaginary / b);
            return c;
        }

        public static Complex operator /(double a, Complex b)
        {
            double denom = b.Real * b.Real + b.Imaginary * b.Imaginary;
            Complex c = new Complex((a * b.Real) / denom, (-a * b.Imaginary) / denom);
            return c;
        }

        public static bool operator ==(Complex a, Complex b)
        {
            return a.Real == b.Real && a.Imaginary == b.Imaginary;
        }

        public static bool operator ==(Complex a, double b)
        {
            return a.Real == b && a.Imaginary == 0;
        }

        public static bool operator ==(double a, Complex b)
        {
            return a == b.Real && 0 == b.Imaginary;
        }

        public static bool operator !=(Complex a, Complex b)
        {
            return a.Real != b.Real || a.Imaginary != b.Imaginary;
        }

        public static bool operator !=(Complex a, double b)
        {
            return a.Real != b || a.Imaginary != 0;
        }

        public static bool operator !=(double a, Complex b)
        {
            return a != b.Real || 0 != b.Imaginary;
        }

        public static implicit operator Complex(double a)
        {
            return new Complex(a, 0);
        }

        public override int GetHashCode()
        {
            return Imaginary == 0? Real.GetHashCode() : Math.Pow(Real, Imaginary).GetHashCode();
        }

        public override bool Equals(object obj)
        {
            switch (obj.GetType().FullName)
            {
                case "CSBoost.Complex" :  return ((Complex)obj) == this;
                case "System.Double" :
                case "System.Single":
                case "System.Byte":
                case "System.SByte":
                case "System.Int16":
                case "System.Int32":
                case "System.Int64":
                case "System.UInt16":
                case "System.UInt32":
                case "System.UInt64":
                case "System.Decimal":
                    return ((double)obj) == this;
            }
            return false;
        }

        public static double Abs(Complex c)
        {
            return XMath.hypot(c.Real, c.Imaginary);
        }

        public static double Arg(Complex c)
        {
            return Math.Atan2(c.Imaginary, c.Real);
        }

        public static Complex Acos(Complex c)
        {
            Complex i1 = new Complex(0, 1);
            Complex r1 = new Complex(1, 0);
            Complex result = -i1 * Log(c + i1 * Sqrt(r1 - c * c));
            return result;
        }

        public static Complex Acosh(Complex c)
        {
            Complex result = Log(c + Sqrt(c + 1) * Sqrt(c - 1));
            return result;
        }

        public static Complex Acot(Complex c)
        {
            Complex i1 = new Complex(0, 1);
            Complex r1 = new Complex(1, 0);
            Complex result = (i1 / new Complex(2, 0)) * (Log(r1 - i1 / c) - Log(r1 + i1 / c));
            return result;
        }

        public static Complex Acoth(Complex c)
        {
            Complex result = 0.5 * (Log(1 + 1/c) - Log(1 - 1/c));
            return result;
        }

        public static Complex Acsc(Complex c)
        {
            Complex i1 = new Complex(0, 1);
            Complex result = -i1 * Log(Sqrt(1 - 1 / (c*c)) + i1/c);
            return result;
        }

        public static Complex Acsch(Complex c)
        {
            Complex result = Log(Sqrt(1 + 1 / (c * c)) + 1 / c);
            return result;
        }

        public static Complex Asec(Complex c)
        {
            Complex i1 = new Complex(0, 1);
            Complex r1 = new Complex(1, 0);
            Complex result = -i1 * Log(Sqrt(r1 / (c * c) - r1) + r1 / c);
            return result;
        }

        public static Complex Asech(Complex c)
        {
            Complex result = Log(1/c + Sqrt(1/c - 1) * Sqrt(1/c + 1));
            return result;
        }

        public static Complex Asin(Complex c)
        {
            Complex i1 = new Complex(0, 1);
            Complex r1 = new Complex(1, 0);
            Complex result = -i1 * Log(c * i1 + Sqrt(r1 - c * c));
            return result;
        }

        public static Complex Asinh(Complex c)
        {
            Complex result = Log(c + Sqrt(c*c + 1));
            return result;
        }

        public static Complex Atan(Complex c)
        {
            Complex i1 = new Complex(0, 1);
            Complex r1 = new Complex(1, 0);
            Complex result = (i1 / new Complex(2, 0)) * (Log(r1 - i1 * c) - Log(r1 + i1 * c));
            return result;
        }

        public static Complex Atanh(Complex c)
        {
            Complex result = 0.5 * (Log(1+c)-Log(1-c));
            return result;
        }

        public static Complex Conjugate(Complex c)
        {
            return new Complex(c.Real, -c.Imaginary);
        }

        public static Complex Cos(Complex c)
        {
            return new Complex(Math.Cos(c.Real) * Math.Cosh(c.Imaginary), -(Math.Sin(c.Real) * Math.Sinh(c.Imaginary)));
        }

        public static Complex Csc(Complex c)
        {
            return 1.0 / Sin(c);
        }

        public static Complex Csch(Complex c)
        {
            return 1.0 / Sinh(c);
        }

        public static Complex Cosh(Complex c)
        {
            return new Complex(Math.Cosh(c.Real) * Math.Cos(c.Imaginary), Math.Sinh(c.Real) * Math.Sin(c.Imaginary));
        }

        public static Complex Cot(Complex c)
        {
            return 1.0 / Tan(c);
        }

        public static Complex Coth(Complex c)
        {
            return 1.0 / Tanh(c);
        }

        public static Complex Exp(Complex c)
        {
            return FromPolarCoordinates(Math.Exp(c.Real), c.Imaginary); 
        }

        public static Complex Log(Complex c)
        {
            return new Complex(Math.Log(c.Magnitude), c.Phase);
        }

        public static Complex Log(Complex a, Complex b)
        {
            return Log(a) / Log(b);
        }

        public static Complex Log(Complex c, double d)
        {
            return Log(c) / Math.Log(d);
        }

        public static Complex Log(double d, Complex c)
        {
            return Math.Log(d) / Log(c);
        }

        public static Complex Log10(Complex c)
        {
            return Log(c) / Math.Log(10);
        }

        public static Complex Pow(Complex a, Complex b)
        {
            return Exp(b * Log(a));
        }

        public static Complex Pow(Complex a, double b)
        {
            return Exp(b * Log(a));
        }

        public static Complex Pow(double a, Complex b)
        {
            return Exp(b * Math.Log(a));
        }

        public static Complex Sec(Complex c)
        {
            return 1.0 / Cos(c);
        }

        public static Complex Sech(Complex c)
        {
            return 1.0 / Cosh(c);
        }

        public static Complex Sin(Complex c)
        {
            return new Complex(Math.Sin(c.Real) * Math.Cosh(c.Imaginary), Math.Cos(c.Real) * Math.Sinh(c.Imaginary));
        }

        public static Complex Sinh(Complex c)
        {
            return new Complex(Math.Sinh(c.Real) * Math.Cos(c.Imaginary), Math.Cosh(c.Real) * Math.Sin(c.Imaginary));
        }

        public static Complex Sqrt(Complex c)
        {
            return Complex.FromPolarCoordinates(Math.Sqrt(c.Magnitude), c.Phase / 2);
        }

        public static Complex Tan(Complex c)
        {
            return Sin(c) / Cos(c);
        }

        public static Complex Tanh(Complex c)
        {
            return Sinh(c) / Cosh(c);
        }

        public override string ToString()
        {
            if (this == 0.0) return "0";
            string result = "", infix = "";
            if (Real != 0)
            {
                result = Real.ToString();
                infix = Imaginary > 0 ? " + " : " - ";
            }
            else if (Imaginary < 0) infix = "-";
            if (Imaginary != 0) result += infix + Math.Abs(Imaginary) + "i";
            return result;
        }
    }
}
