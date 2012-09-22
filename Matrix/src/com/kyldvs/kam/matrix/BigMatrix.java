package com.kyldvs.kam.matrix;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import com.kyldvs.kam.matrix.misc.InvalidDimensionsException;

public class BigMatrix {

	public static BigMatrix create(int[][] arr) {
		BigInteger[][] m = new BigInteger[arr.length][arr[0].length];
		for (int r = 0; r < arr.length; r++) {
			for (int c = 0; c < arr[0].length; c++) {
				m[r][c] = BigInteger.valueOf(arr[r][c]);
			}
		}
		return new BigMatrix(m, true);
	}
	
	public static BigMatrix create(long[][] arr) {
		BigInteger[][] m = new BigInteger[arr.length][arr[0].length];
		for (int r = 0; r < arr.length; r++) {
			for (int c = 0; c < arr[0].length; c++) {
				m[r][c] = BigInteger.valueOf(arr[r][c]);
			}
		}
		return new BigMatrix(m, true);
	}
	
	public static BigMatrix power(BigMatrix m, int pow) {
		if (pow < 0) {
			throw new RuntimeException("Can not use a negative exponent on integers.");
		}
		if (m.cols() != m.rows()) {
			throw new InvalidDimensionsException();
		}
		return power(m, pow, 1, new HashMap<Integer, BigMatrix>());
	}
	
	private static BigMatrix power(BigMatrix m, int pow, int at, Map<Integer, BigMatrix> map) {
		if (map.containsKey(at * pow)) {
			return map.get(at * pow);
		}
		
		BigMatrix ret;
		if (pow == 0) {
			int[][] arr = new int[m.rows()][m.cols()];
			for (int i = 0; i < arr.length; i++) {
				Arrays.fill(arr[i], 1);
			}
			ret = BigMatrix.create(arr);
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
	
	private static BigInteger[][] copy(BigInteger[][] arr) {
		BigInteger[][] ret = new BigInteger[arr.length][arr[0].length];
		for (int i = 0; i < arr.length; i++) {
			ret[i] = Arrays.copyOf(arr[i], arr[i].length);
		}
		return ret;
	}
	
	protected BigInteger[][] m;

	public BigMatrix() {
		m = new BigInteger[0][0];
	}
	
	public BigMatrix(BigInteger[][] m) {
		this(m, false);
	}
	
	private BigMatrix(BigInteger[][] m, boolean fast) {
		if (fast) {
			this.m = m;
		} else {
			this.m = copy(m);
		}
	}
	
	public final BigInteger get(int r, int c) {
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
	
	public BigMatrix multiply(BigMatrix m) {
		if (cols() != m.rows()) {
			throw new InvalidDimensionsException();
		}
		
		BigInteger[][] result = new BigInteger[rows()][m.cols()];
		for (int r = 0; r < rows(); r++) {
			for (int c = 0; c < m.cols(); c++) {
				result[r][c] = BigInteger.ZERO;
				for (int i = 0; i < cols(); i++) {
					result[r][c] = result[r][c].add(get(r, i).multiply(m.get(i, c)));
				}
			}
		}

		return new BigMatrix(result, true);
	}
	
	@Override
	public boolean equals(Object o) {
		if (o instanceof BigMatrix) {
			return Arrays.deepEquals(m, ((BigMatrix) o).m);
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
				strs[r][c] = m[r][c].toString();
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
				if (c + 1 < cols()) {
					sb.append(", ");
				}
			}
			sb.append("]\n");
		}
		return sb.toString();
	}
	
	public static void main(String[] args) {
		BigMatrix fib3 = BigMatrix.create(new int[][]{
			{1,1,1},
			{-3,0,0},
			{0,2,0}
		});

		BigMatrix start = BigMatrix.create(new int[][]{
				{0,0,1},
		});
		
		int n = 100;
		for (int i = 0; i < n; i++) {
			BigMatrix res = start.multiply(BigMatrix.power(fib3, i));
			System.out.println("F_" + i + ": " + res.get(0, 0));
		}
	}
}
