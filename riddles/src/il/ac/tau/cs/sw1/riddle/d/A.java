package il.ac.tau.cs.sw1.riddle.d;

public class A {

	private int l = 7;

	public static void main(String[] args) {
		B b1 = new B();//2
		B b2 = new B();//3
		B b3 = new B();//7
		B b4 = new B();//5
		setTheBs(b1, b2, b3, b4);
		System.out.println(b1.getI() * b2.getI() * b3.getI()
				* b4.getI());
	}

	private static void setTheBs(B b1, B b2, B b3, B b4) {
		int k = 5;
		b1.setI(B.I); // replace
		for (int j = 0; j < 11; j++) {
			if (j > b1.getI()
					&& (j == 2 || j == 3 || j == 5 || j == 7)) {
				b2.setI(j); // replace
				break;
			}
		}//j==7
		A a = new A();
		b3.setI(a.l); // replace
		a.instanceSetTheBs(b1, b3, b4, k);
	}

	private void instanceSetTheBs(B b1, B b3, B b4, int k) {
		if (b3.getI() > b1.getI()) {
			b4.setI(k); // replace
		}
	}
}