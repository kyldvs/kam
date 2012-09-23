package com.kyldvs.kam.matrix;

import static org.junit.Assert.*;

import org.junit.Test;

import com.kyldvs.kam.matrix.misc.InvalidDimensionsException;

public class MatrixTests {

	private static final Matrix a = new Matrix(new double[][] {
			{1, 1},
			{1, 0}
	});
	
	private static final Matrix b = new Matrix(new double[][] {
			{1, 1, 2},
			{1, 0, 2}
	});
	
	
	@Test
	public void power1() {
		Matrix act = Matrix.power(a, 10);
		Matrix exp = new Matrix(new double[][] {
				{89, 55},
				{55, 34}
		});
		assertEquals(exp, act);
	}
	
	@Test(expected = InvalidDimensionsException.class)
	public void power2() {
		Matrix.power(Matrix.zero(5, 6), 10);
	}
	
	@Test
	public void solve() {
		Matrix act = Matrix.solve(b);
		Matrix exp = new Matrix(new double[][] {
				{1, 0, 2},
				{0, 1, 0}
		});
		assertEquals(exp, act);
	}

}
