package com.kyldvs.kam.matrix;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import com.kyldvs.kam.matrix.misc.InvalidDimensionsException;

public class Matrix {

	public static int DECIMALS = 2;
	
	public static Matrix power(Matrix m, int pow) {
		if (m.cols() != m.rows()) {
			throw new InvalidDimensionsException();
		}
		if (pow < 0) {
			Matrix pos = power(m, -pow, 1, new HashMap<Integer, Matrix>());
			double[][] neg = new double[pos.rows()][pos.cols()];
			for (int r = 0; r < neg.length; r++) {
				for (int c = 0; c < neg[0].length; c++) {
					neg[r][c] = 1.0 / pos.get(r, c);
				}
			}
			return new Matrix(neg, true);
		}
		return power(m, pow, 1, new HashMap<Integer, Matrix>());
	}
	
	private static Matrix power(Matrix m, int pow, int at, Map<Integer, Matrix> map) {
		if (map.containsKey(pow)) {
			return map.get(pow);
		}
		
		Matrix ret;
		if (pow == 0) {
			double[][] arr = new double[m.rows()][m.cols()];
			for (int i = 0; i < arr.length; i++) {
				Arrays.fill(arr[i], 1D);
			}
			ret = new Matrix(arr);
		} else if (pow == 1) {
			ret = m;
		} else if (pow == 2) {
			ret = m.multiply(m);
		} else if (pow % 2 == 0) {
			ret = power(power(m, pow/2, at, map), 2, pow/2, map);
		} else {
			ret = power(power(m, pow/2, at, map), 2, pow/2, map).multiply(m);
		}
		
		map.put(at * pow, ret);
		return ret;
	}
	
	private static double[][] copy(double[][] arr) {
		double[][] ret = new double[arr.length][arr[0].length];
		for (int i = 0; i < arr.length; i++) {
			ret[i] = Arrays.copyOf(arr[i], arr[i].length);
		}
		return ret;
	}

	protected double[][] m;

	public Matrix() {
		m = new double[0][0];
	}
	
	public Matrix(double[][] m) {
		this(m, false);
	}
	
	private Matrix(double[][] m, boolean fast) {
		if (fast) {
			this.m = m;
		} else {
			this.m = copy(m);
		}
	}
	
	public final double get(int r, int c) {
		if (r < 0 || r >= rows() || c < 0 || c >= cols()) {
			throw new IndexOutOfBoundsException();
		}
		return m[r][c];
	}
	
	public final int rows() {
		return m.length;
	}
	
	public final int cols() {
		return m[0].length;
	}
	
	public Matrix multiply(Matrix m) {
		if (cols() != m.rows()) {
			throw new InvalidDimensionsException();
		}
		
		double[][] result = new double[rows()][m.cols()];
		for (int r = 0; r < rows(); r++) {
			for (int c = 0; c < m.cols(); c++) {
				for (int i = 0; i < cols(); i++) {
					result[r][c] += get(r, i) * m.get(i, c);
				}
			}
		}

		return new Matrix(result, true);
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof Matrix) {
			return Arrays.deepEquals(m, ((Matrix) o).m);
		}
		return false;
	}
	
	@Override
	public int hashCode() {
		return Arrays.deepHashCode(m);
	}
	
	@Override
	public String toString() {
		String[][] strs = new String[rows()][cols()];
		int[] len = new int[cols()];
		for (int r = 0; r < rows(); r++) {
			for (int c = 0; c < cols(); c++) {
				strs[r][c] = String.format("%." + DECIMALS + "f", m[r][c]);
				len[c] = Math.max(len[c], strs[r][c].length());
			}
		}
		
		StringBuilder sb = new StringBuilder();
		for (int r = 0; r < rows(); r++) {
			sb.append("[");
			for (int c = 0; c < cols(); c++) {
				for (int ct = 0; ct < len[c] - strs[r][c].length(); ct++) {
					sb.append(" ");
				}
				sb.append(strs[r][c]);
				if (c + 1 < rows()) {
					sb.append(", ");
				}
			}
			sb.append("]\n");
		}
		return sb.toString();
	}
	
	public static void main(String[] args) {
		Matrix.DECIMALS = 1;
		
		Matrix fib = new Matrix(new double[][]{
				{1,1},
				{1,0}
			});
		
		Matrix res = Matrix.power(fib, 100);
		System.out.println(res.toString());
	}
}
