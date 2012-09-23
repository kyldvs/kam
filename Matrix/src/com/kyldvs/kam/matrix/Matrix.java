package com.kyldvs.kam.matrix;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import com.google.common.base.Function;
import com.kyldvs.kam.matrix.misc.InvalidDimensionsException;

/**
 * An immutable matrix class, with utilities
 * 
 * @author kyle
 */
public class Matrix {

	private static final double EPS = 1e-11;

	/**
	 * A function to normalize everything that should be treated as 0
	 */
	public static final Function<Double, Double> ZERO = new Function<Double, Double>() {
		@Override
		public Double apply(Double obj) {
			return Math.abs(obj) < EPS ? 0 : obj;
		}
	};

	/**
	 * @param m
	 * @param f
	 * @return applies f to each entry of m and returns the resulting matrix
	 */
	public static Matrix apply(Matrix m, Function<Double, Double> f) {
		double[][] arr = new double[m.rows()][m.cols()];
		for (int row = 0; row < m.rows(); row++) {
			for (int col = 0; col < m.cols(); col++) {
				arr[row][col] = f.apply(m.fastGet(row, col));
			}
		}
		return new Matrix(arr, true);
	}
	
	public static Matrix create(int r, int c, double d) {
		double[][] arr = new double[r][c];
		for (int i = 0; i < r; i++) {
			Arrays.fill(arr[i], d);
		}
		return new Matrix(arr, true);
	}

	public static Matrix create(Matrix m) {
		return new Matrix(m.arr, true);
	}

	public static Matrix zero(int r, int c) {
		return new Matrix(new double[r][c], true);
	}

	public static Matrix one(int r, int c) {
		return create(r, c, 1.0);
	}

	public static Matrix identity(int n) {
		double[][] arr = new double[n][n];
		for (int i = 0; i < n; i++) {
			arr[i][i] = 1.0;
		}
		return new Matrix(arr, true);
	}

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

	private static Matrix power(Matrix m, int pow, int at,
			Map<Integer, Matrix> map) {
		if (map.containsKey(at * pow)) {
			return map.get(at * pow);
		}

		Matrix ret;
		if (pow == 0) {
			ret = one(m.rows(), m.cols());
		} else if (pow == 1) {
			ret = m;
		} else if (pow == 2) {
			ret = m.multiply(m);
		} else if (pow % 2 == 0) {
			ret = power(power(m, pow / 2, at, map), 2, pow / 2, map);
		} else {
			ret = power(power(m, pow / 2, at, map), 2, pow / 2, map)
					.multiply(m);
		}

		map.put(at * pow, ret);
		return ret;
	}

	public static double determinant(Matrix m) {
		if (m.rows() != m.cols()) {
			throw new InvalidDimensionsException();
		}

		double[][] identity = identity(m.rows()).arr;
		double[][] arr = copy(m.arr);
		double[][] upperTriangular = upperTriangular2(arr, identity);

		int swaps = 0;
		double prod = 1;
		for (int i = 0; i < upperTriangular.length; i++) {
			if (Math.abs(identity[i][i]) < EPS) {
				swaps++;
			}
			prod *= upperTriangular[i][i];
		}

		if (swaps % 2 == 1) {
			prod = -prod;
		}

		return prod;
	}

	public static Matrix inverse(Matrix m) {
		if (m.rows() != m.cols()) {
			throw new InvalidDimensionsException();
		}

		double[][] decomp = identity(m.rows()).arr;
		double[][] arr = copy(m.arr);
		double[][] upperTriangular = upperTriangular2(arr, decomp);
		double[][] specialUpperTriangular = special2(upperTriangular, decomp);
		backSolve2(specialUpperTriangular, decomp);
		return new Matrix(decomp);
	}

	/**
	 * @param m
	 * @return the permuted lower triangular matrix in the LU decomposition of m
	 */
	public static Matrix LUDecomposition(Matrix m) {
		double[][] decomp = identity(m.rows()).arr;
		double[][] arr = copy(m.arr);
		upperTriangular2(arr, decomp);
		return new Matrix(decomp);
	}

	public static Matrix solve(Matrix m) {
		double[][] arr = copy(m.arr);
		double[][] upperTriangular = upperTriangular(arr);
		double[][] specialUpperTriangular = special(upperTriangular);
		double[][] backSolve = backSolve(specialUpperTriangular);

		return new Matrix(backSolve, true);
	}

	public static Matrix upperTriangular(Matrix m) {
		double[][] arr = copy(m.arr);
		return new Matrix(upperTriangular(arr), true);
	}

	public static Matrix specialUpperTriangular(Matrix m) {
		double[][] arr = copy(m.arr);
		return new Matrix(special(upperTriangular(arr)), true);
	}

	private static double[][] upperTriangular(double[][] arr) {
		return upperTriangular2(arr, null);
	}

	private static double[][] upperTriangular2(double[][] arr, double[][] copy) {
		boolean validCopy = copy != null && arr.length == copy.length;

		int rows = arr.length;
		int cols = arr[0].length;

		int row = 0;
		for (int col = 0; col < cols; col++) {
			if (row >= rows) {
				break;
			}

			int pivot = row;
			double pivotValue = Math.abs(arr[row][col]);
			for (int i = row; i < rows; i++) {
				double value = Math.abs(arr[i][col]);
				if (value < EPS) {
					arr[i][col] = 0;
				} else if (value > pivotValue) {
					pivot = i;
					pivotValue = value;
				}
			}

			if (pivotValue < EPS) {
				continue;
			}

			swapRows(row, pivot, arr);
			if (validCopy) {
				swapRows(row, pivot, copy);
			}

			for (int j = row + 1; j < rows; j++) {
				double coeff = -arr[j][col] / arr[row][col];

				addRows(row, coeff, j, arr);
				if (validCopy) {
					addRows(row, coeff, j, copy);
				}
			}

			row++;
		}
		return arr;
	}

	private static double[][] special(double[][] arr) {
		return special2(arr, null);
	}

	private static double[][] special2(double[][] arr, double[][] copy) {
		boolean validCopy = copy != null && arr.length == copy.length;

		int rows = arr.length;
		int cols = arr[0].length;

		for (int row = 0; row < rows; row++) {
			int nonZero = -1;
			for (int col = 0; col < cols; col++) {
				if (Math.abs(arr[row][col]) > EPS) {
					nonZero = col;
					break;
				}
			}

			if (nonZero == -1) {
				continue;
			}

			double coeff = 1.0 / arr[row][nonZero];
			multiplyRow(row, coeff, arr);
			if (validCopy) {
				multiplyRow(row, coeff, copy);
			}
		}
		return arr;
	}

	private static double[][] backSolve(double[][] arr) {
		return backSolve2(arr, null);
	}

	private static double[][] backSolve2(double[][] arr, double[][] copy) {
		boolean validCopy = copy != null && arr.length == copy.length;

		int rows = arr.length;
		int cols = arr[0].length;

		for (int row = rows - 1; row >= 0; row--) {
			int nonZero = -1;
			for (int col = 0; col < cols; col++) {
				if (Math.abs(arr[row][col]) > EPS) {
					nonZero = col;
					break;
				}
			}

			if (nonZero == -1) {
				continue;
			}

			for (int r2 = row - 1; r2 >= 0; r2--) {
				double coeff = -arr[r2][nonZero] / arr[row][nonZero];
				addRows(row, coeff, r2, arr);
				if (validCopy) {
					addRows(row, coeff, r2, copy);
				}
			}
		}

		return arr;
	}

	private static void swapRows(int r1, int r2, double[][]... arrs) {
		for (double[][] arr : arrs) {
			double[] tmp = arr[r1];
			arr[r1] = arr[r2];
			arr[r2] = tmp;
		}
	}

	private static void addRows(int row, double coeff, int to,
			double[][]... arrs) {
		for (double[][] arr : arrs) {
			for (int col = 0; col < arr[row].length; col++) {
				arr[to][col] += coeff * arr[row][col];
			}
		}
	}

	private static void multiplyRow(int row, double coeff, double[][]... arrs) {
		for (double[][] arr : arrs) {
			for (int col = 0; col < arr[row].length; col++) {
				arr[row][col] *= coeff;
			}
		}
	}

	private static double[][] copy(double[][] arr) {
		double[][] ret = new double[arr.length][arr[0].length];
		for (int i = 0; i < arr.length; i++) {
			ret[i] = Arrays.copyOf(arr[i], arr[i].length);
		}
		return ret;
	}

	/*
	 * Matrix implementation below this.
	 */

	private double[][] arr;

	public Matrix() {
		arr = new double[0][0];
	}

	public Matrix(double[][] arr) {
		this(arr, false);
	}

	/**
	 * Private constructor that gives the option of referencing the array
	 * directly, rather than copying its contents
	 * 
	 * @param arr
	 * @param fast
	 */
	private Matrix(double[][] arr, boolean fast) {
		if (fast) {
			this.arr = arr;
		} else {
			this.arr = copy(arr);
		}
	}

	public final double get(int r, int c) {
		if (r < 0 || r >= rows() || c < 0 || c >= cols()) {
			throw new IndexOutOfBoundsException();
		}
		return arr[r][c];
	}

	private final double fastGet(int r, int c) {
		return arr[r][c];
	}

	public final int rows() {
		return arr.length;
	}

	public final int cols() {
		return arr[0].length;
	}

	public Matrix add(Matrix m) {
		if (cols() != m.cols() || rows() != m.rows()) {
			throw new InvalidDimensionsException();
		}

		double[][] result = new double[rows()][cols()];
		for (int r = 0; r < rows(); r++) {
			for (int c = 0; c < cols(); c++) {
				result[r][c] = fastGet(r, c) + m.fastGet(r, c);
			}
		}

		return new Matrix(result, true);
	}

	public Matrix multiply(Matrix m) {
		if (cols() != m.rows()) {
			throw new InvalidDimensionsException();
		}

		double[][] result = new double[rows()][m.cols()];
		for (int r = 0; r < rows(); r++) {
			for (int c = 0; c < m.cols(); c++) {
				for (int i = 0; i < cols(); i++) {
					result[r][c] += fastGet(r, i) * m.fastGet(i, c);
				}
			}
		}

		return new Matrix(result, true);
	}

	private Matrix transpose = null;

	public Matrix transpose() {
		if (transpose != null) {
			return transpose;
		}

		double[][] result = new double[cols()][rows()];
		for (int row = 0; row < rows(); row++) {
			for (int col = 0; col < cols(); col++) {
				result[col][row] = fastGet(row, col);
			}
		}

		transpose = new Matrix(result, true);
		return transpose;
	}

	private Boolean symmetric = null;

	public boolean isSymmetric() {
		if (symmetric != null) {
			return symmetric;
		}
		symmetric = this.equals(transpose());
		return symmetric;
	}

	private Integer rank = null;

	public int rank() {
		if (rank != null) {
			return rank;
		}

		Matrix m = upperTriangular(this);
		int pivots = 0;
		for (int row = 0; row < m.rows(); row++) {
			for (int col = 0; col < m.cols(); col++) {
				if (Math.abs(m.fastGet(row, col)) > EPS) {
					pivots++;
					break;
				}
			}
		}

		rank = pivots;
		return rank;
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof Matrix) {
			double[][] tht = ((Matrix) o).arr;
			if (arr.length != tht.length || arr[0].length != tht[0].length) {
				return false;
			}
			for (int i = 0; i < arr.length; i++) {
				for (int j = 0; j < arr[0].length; j++) {
					if (Math.abs(arr[i][j] - tht[i][j]) > EPS) {
						return false;
					}
				}
			}
			return true;
		}
		return false;
	}

	@Override
	public int hashCode() {
		return Arrays.deepHashCode(arr);
	}

	public static int DECIMALS = 2;

	@Override
	public String toString() {
		int decimals = Math.max(0, DECIMALS);
		String[][] strs = new String[rows()][cols()];
		int[] len = new int[cols()];
		for (int r = 0; r < rows(); r++) {
			for (int c = 0; c < cols(); c++) {
				strs[r][c] = String
						.format("%." + decimals + "f", fastGet(r, c));
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
			sb.append("]");
			if (r + 1 < rows()) {
				sb.append("\n");
			}
		}
		return sb.toString();
	}

	/*
	 * Simple main method
	 */

	public static void main(String[] args) {
		Matrix.DECIMALS = 2;

		Matrix a = new Matrix(new double[][] { { 1, -1, 7, -7 },
				{ 1, 2, 3, -3 }, { 2, 3, 4, 3 } });
		System.out.println("Matrix A:");
		System.out.println(a);
		System.out.println("\nA Transpose:");
		System.out.println(a.transpose());
		System.out.println("\nA is Symmetric: " + a.isSymmetric());
		System.out.println("\nSolution to A:");
		System.out.println(solve(a));
		Matrix LUD = LUDecomposition(a);
		System.out.println("\nLUDecomposition of A:");
		System.out.println(LUD);
		System.out.println("\nLUD.A:");
		System.out.println(LUD.multiply(a));
		System.out.println("\nUpper Triangular A:");
		System.out.println(upperTriangular(a));

		Matrix b = new Matrix(new double[][] { { 1, 2 }, { 0, 1 } });
		System.out.println("\nMatrix B:");
		System.out.println(b);
		Matrix bi = inverse(b);
		System.out.println("\nB Inverse:");
		System.out.println(bi);
		System.out.println("\nB.BInverse:");
		System.out.println(b.multiply(bi));
		System.out.println("\nBInverse.B:");
		System.out.println(bi.multiply(b));

		Matrix c = new Matrix(new double[][] { { 1, 3, 5 }, { 3, -5, 7 },
				{ 5, 7, 2 } });

		System.out.println("\nMatrix C:");
		System.out.println(c);
		System.out.println("\nC Transpose:");
		System.out.println(c.transpose());
		System.out.println("\nC is symmetric: " + c.isSymmetric());

		Matrix d = new Matrix(new double[][] { { 1, 3, 2, -1, 0 },
				{ 2, 6, 1, 4, 3 }, { -1, -3, -3, 3, 1 }, { 3, 9, 8, -7, 2 } });

		System.out.println("\nMatrix D:");
		System.out.println(d);
		System.out.println("\nUpper Triangular D:");
		System.out.println(upperTriangular(d));
		System.out.println("\nD Rank: " + d.rank());
	}
}
