package il.ac.tau.cs.sw1.hw6;

public class Polynomial
{//[a,b,c,..,z] == a + b*x + c*x^2 + ... + z*x^21
	private double[] polyArr;
	/*
	 * Creates the zero-polynomial with p(x) = 0 for all x.
	 */
	public Polynomial()
	{
		this.polyArr = new double[1];
		this.polyArr[0] = 0.0; //default is already 0.0, added for readability
	} 
	/*
	 * Creates a new polynomial with the given coefficients.
	 */
	public Polynomial(double[] coefficients) 
	{
		this.polyArr = new double[coefficients.length];
		for(int i=0; i<polyArr.length; i++)
		{
			this.polyArr[i] = coefficients[i];
		}
	}
	/*
	 * Addes this polynomial to the given one
	 *  and retruns the sum as a new polynomial.
	 */
	public Polynomial adds(Polynomial polynomial)
	{
		if(this.polyArr.length >= polynomial.polyArr.length)
		{
			Polynomial retPoly = new Polynomial(this.polyArr);
			return addLoop(retPoly, polynomial);
		}
		else
		{
			Polynomial retPoly = new Polynomial(polynomial.polyArr);
			return addLoop(retPoly, this);
		}		
	}
	public Polynomial addLoop(Polynomial retPoly, Polynomial shortPoly)
	{
		for(int i=0; i<shortPoly.polyArr.length; i++)
		{
			retPoly.polyArr[i] += shortPoly.polyArr[i];
		}
		return retPoly;	
	}
	
	/*
	 * Multiplies a to this polynomial and returns 
	 * the result as a new polynomial.
	 */
	public Polynomial multiply(double a)
	{
		double[] coef = new double[this.polyArr.length];
		for(int i=0; i<coef.length; i++)
		{
			coef[i] = a * this.polyArr[i];
			if(coef[i] == -0.0)
				coef[i] = 0.0;
		}
		return new Polynomial(coef);
	}
	/*
	 * Returns the degree (the largest exponent) of this polynomial.
	 */
	public int getDegree()
	{
		for(int i=this.polyArr.length-1; i>=0; i--)
		{
			if(this.polyArr[i] != 0)
			{
				return i;
			}
		}
		return 0;
	}
	/*
	 * Returns the coefficient of the variable x 
	 * with degree n in this polynomial.
	 */
	public double getCoefficient(int n)
	{
		if(n >= this.polyArr.length)
		{
			return 0.0;
		}
		else
		{
			return this.polyArr[n];
		}
	}
	
	/*
	 * set the coefficient of the variable x 
	 * with degree n to c in this polynomial.
	 * If the degree of this polynomial < n, it means that that the coefficient of the variable x 
	 * with degree n was 0, and now it will change to c. 
	 */
	public void setCoefficient(int n, double c)
	{
		if(n > this.getDegree()) 
		{
			double[] coef = new double[n+1];
			for (int i=0; i<this.polyArr.length; i++) 
			{
				coef[i] = this.polyArr[i];
			}
			this.polyArr = coef;	
		}
		this.polyArr[n] = c;
	}
	
	/*
	 * Returns the first derivation of this polynomial.
	 *  The first derivation of a polynomal a0x0 + ...  + anxn is defined as 1 * a1x0 + ... + n anxn-1.
	
	 */
	public Polynomial getFirstDerivation()
	{
		if(this.getDegree() == 0)
		{
			return new Polynomial();
		}
		else
		{
			double[] deriv = new double[this.getDegree()];
			for(int i=0; i<deriv.length; i++)
			{
				deriv[i] = (i+1) * this.polyArr[i+1];
			}
			return new Polynomial(deriv);
		}
	}
	
	/*
	 * given an assignment for the variable x,
	 * compute the polynomial value
	 */
	public double computePolynomial(double x)
	{
		double val = this.polyArr[0];
		for(int i=1; i<this.polyArr.length; i++)
		{
			val += this.polyArr[i] * Math.pow(x, i);
		}
		return val;
	}
	
	/*
	 * given an assignment for the variable x,
	 * return true iff x is an extrema point (local minimum or local maximum of this polynomial)
	 * x is an extrema point if and only if The value of first derivation of a polynomal at x is 0
	 * and the second derivation of a polynomal value at x is not 0.
	 */
	public boolean isExtrema(double x)
	{
		Polynomial deriv1 = this.getFirstDerivation();
		Polynomial deriv2 = deriv1.getFirstDerivation();
		return deriv1.computePolynomial(x) == 0.0 && deriv2.computePolynomial(x) != 0.0;
	}
	
	
	
	

    
    

}
