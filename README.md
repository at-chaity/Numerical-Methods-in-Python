# Numerical-Methods-in-Python

This repository presents a comprehensive and modular implementation of **Numerical Methods**, covering key topics in scientific computing. Each method is implemented as a Python function with accompanying **Jupyter notebooks** that demonstrate accuracy, efficiency, and visual comparisons using real problems.

This work reflects both algorithmic understanding and practical numerical analysis skills, developed through structured coding, symbolic computation, and clear presentation.

The codes were built following the algorithms of book **Numerical Analysis (9th Edition) -Richard L. Burden, J. Douglas Faires**.

---

## 🧠 Topics Covered

### 1. Root Finding
- **Methods:** Bisection, Newton-Raphson
- **Features:** Symbolic parsing, convergence control, accuracy tracking
- **Notebook:** Side-by-side comparisons

### 2. Numerical Differentiation
- **Methods:** Forward, Backward, Three Point Mid Point, Three Point End Point , Five-Point End Point
- **Notebook:** Method-wise accuracy tested on data

### 3. Interpolation
- **Methods:** Lagrange, Newton's Divided Differences, Cubic Spline (Manual & `scipy`)
- **Notebook:** Function reconstruction and comparative plots

### 4. Numerical Integration
- **Methods:** Composite Simpson’s Rule, Trapezoidal Rule, Romberg Integration
- **Extras:** Double integral using Simpson’s rule vs `scipy.integrate.dblquad`
- **Notebook:** Tabular and visual error analysis

### 5. Initial Value Problems (IVP)
- **Methods:** Euler's Method, Runge-Kutta order 4, Predictor-Corrector
- **Notebook:** Side-by-side performance and convergence test

### 6. Boundary Value Problems (BVP)
- **Method:** Linear Shooting Method
- **Notebook:** Comparison with `scipy.integrate.solve_bvp`

### 7. Linear Systems
- **Methods:** Gaussian Elimination, LU Decomposition, Gauss-Seidel, Jacobi, SOR, SVD
- **Notebook:** Accuracy, iteration count, and stability comparison

### 8. Nonlinear Systems
- **Methods:** Newton’s Method, Steepest Descent
- **Notebook:** Convergence behavior and error tracking

---

## 📊 Features

- **Symbolic Function Parsing:** Using `sympy` to handle user-defined functions
- **Tabulated Results:** Using `pandas` for clean formatting
- **Visual Insights:** Error plots and convergence curves with `matplotlib`
- **Clean Modularity:** All core methods are written as reusable, testable Python functions

---

## 🛠 Libraries Used

- `numpy`
- `scipy`
- `sympy`
- `matplotlib`
- `pandas`

---

## 🧪 How to Use

1. Clone the repository:
   ```bash
   git clone https://github.com/at-chaity/Numerical-Methods-in-Python.git
   cd numerical-methods
2. Launch Jupyter Notebook.
3. Navigate to the Notebooks/ folder to explore comparisons and visualizations.
---

## 📁 Project Structure

Numerical-Methods/
- Codes/               (All core method implementations)
- Notebooks/           (Interactive comparisons and analysis)
- README.md            (Project documentation)

---

## ✨ Objective
This project was created to solidify core numerical algorithms and showcase the comparisons using Python. It aims to strike a balance between clarity, mathematical rigor, and readability — making it a valuable resource for academic use, self-learning, or teaching assistance.

Crafted with a lot of care and curiosity.
— Afifa Tabassum Chaity


