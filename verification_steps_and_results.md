# 詳細驗證步驟過程及結果 (Detailed Verification Steps and Results)

## 概述 (Overview)

本文檔詳細記錄了 R 語言實現的 `thd2arma` 函數的驗證過程，確保其邏輯與 MATLAB 版本完全匹配。驗證涵蓋從簡單的單變量模型到複雜的多變量季節模型。

## 驗證目標 (Verification Objectives)

- 驗證 R 實現的邏輯與 MATLAB 版本相符
- 確保輸出矩陣的維度正確
- 驗證多項式乘法算法的正確性
- 測試簡單和複雜案例的結果

## 驗證步驟 (Verification Steps)

### 1. 簡單單變量 ARMA(1,1) 模型測試 (Simple Univariate ARMA(1,1) Test)

#### 測試設定
```R
# 簡單的單變量 ARMA(1,1) 模型 - 無季節成分
# Simple univariate ARMA(1,1) model - no seasonal components
m <- 1  # 一個內生變量 / one endogenous variable
r <- 0  # 無外生變量 / no exogenous variables
s <- 4  # 季節週期 / seasonal period (not used in this simple case)
p <- 1  # 常規 AR 階數 / regular AR order
P <- 0  # 季節 AR 階數 / seasonal AR order
q <- 1  # 常規 MA 階數 / regular MA order
Q <- 0  # 季節 MA 階數 / seasonal MA order
g <- 0  # 外生變量滯後階數 / exogenous variable lag order

# 參數設定 / Parameter settings
theta <- c(
  # FR1 (1 參數) - AR 係數 / AR coefficient
  0.5,        # AR coefficient
  # FS (0 參數) - 無季節 AR / no seasonal AR
  # AR1 (1 參數) - MA 係數 / MA coefficient
  0.3,        # MA coefficient
  # G (0 參數) - 無外生 / no exogenous
  # V (1 參數) - 方差 / variance
  1.0         # variance
)
```

#### 期望結果
對於簡單的 ARMA(1,1) 模型 `F(B)y(t) = A(B)e(t)`:
- F(B) = I + F1*B = [1] + [0.5]*B → F = [1, 0.5]
- A(B) = I + A1*B = [1] + [0.3]*B → A = [1, 0.3]

#### 實際執行結果
```
Simple ARMA(1,1) Results:
F Matrix (should be [1, 0.5]):
     [,1] [,2]
[1,]    1  0.5
A Matrix (should be [1, 0.3]):
     [,1] [,2]
[1,]    1  0.3
G Matrix (should be [0]):
    
[1,]
V Matrix (should be [1]):
     [,1]
[1,]    1
```

#### 驗證結果
```
Simple test validation:
F matrix matches expected [1, 0.5]: TRUE 
A matrix matches expected [1, 0.3]: TRUE 
PASS: Simple ARMA(1,1) test matches expected results
```

### 2. 複雜多變量季節模型測試 (Complex Multivariate Seasonal Model Test)

#### 測試設定
```R
# 複雜模型設定
m <- 2  # 兩個內生變量 / two endogenous variables
r <- 1  # 一個外生變量 / one exogenous variable
s <- 4  # 季節性 (季度) / quarterly seasonality
p <- 1  # 常規 AR 階數 / regular AR order
P <- 0  # 季節 AR 階數 / seasonal AR order
q <- 1  # 常規 MA 階數 / regular MA order
Q <- 1  # 季節 MA 階數 / seasonal MA order
g <- 0  # 外生變量滯後階數 / exogenous variable lag order

k <- max(p + P * s, q + Q * s, g)  # max(1, 1+4, 0) = 5
n_states <- k * m  # 10

# 模型參數 / Model parameters
theta <- c(
  # FR1 (4 參數) - AR 係數 / AR coefficients
  0.5, 0.1,  # FR1 第一行 / First row of FR1
  0.2, 0.6,  # FR1 第二行 / Second row of FR1
  # FS (0 參數) - 無季節 AR / Seasonal AR (empty for this test)
  # AR1 (4 參數) - MA 係數 / MA coefficients
  0.3, -0.1, # AR1 第一行 / First row of AR1
  0.1, 0.4,  # AR1 第二行 / Second row of AR1
  # AS1 (4 參數) - 季節 MA 係數 / Seasonal MA coefficients
  0.0, 0.2,  # AS1 第一行 / First row of AS1
  0.2, 0.0,  # AS1 第二行 / Second row of AS1
  # G0 (2 參數)
  1.1, 1.2,
  # V (3 參數，對稱共變異數) / V (3 params for symmetric covariance)
  1.0,   # V[1,1]
  0.5,   # V[2,1] 
  2.0    # V[2,2]
)
```

#### 分解函數結果 (Factored Components)
```
Factored components:
FR: 2 2 
     [,1] [,2]
[1,]  0.5  0.2
[2,]  0.1  0.6
FS: 2 0  (should be 2x0)
    
[1,]
[2,]
AR: 2 2 
     [,1] [,2]
[1,]  0.3  0.1
[2,] -0.1  0.4
AS: 2 2 
     [,1] [,2]
[1,]  0.0  0.2
[2,]  0.2  0.0
```

#### F 矩陣驗證 (F Matrix Verification)
對於 F(B) = FR(B) * FS(B^s) = (I + FR1*B) * (I + 0) = I + FR1*B
因此 F 應該是 [I, FR1] 格式

```
Expected F based on algorithm [I, FR1]:
     [,1] [,2] [,3] [,4]
[1,]    1    0  0.5  0.2
[2,]    0    1  0.1  0.6
Actual F:
     [,1] [,2] [,3] [,4]
[1,]    1    0  0.5  0.2
[2,]    0    1  0.1  0.6
F matrix matches expected [I, FR1]: TRUE 
```

#### A 矩陣驗證 (A Matrix Verification)
對於 A(B) = AR(B) * AS(B^s) = (I + AR1*B) * (I + AS1*B^4)
展開後: I + AR1*B + AS1*B^4 + AR1*AS1*B^5
因此 A 應該包含: [I, AR1, 0, 0, AS1, AR1%*%AS1]

```
Expected A based on polynomial multiplication:
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
[1,]    1    0  0.3  0.1    0    0    0    0  0.0   0.2  0.02  0.06
[2,]    0    1 -0.1  0.4    0    0    0    0  0.2   0.0  0.08 -0.02
Actual A:
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
[1,]    1    0  0.3  0.1    0    0    0    0  0.0   0.2  0.02  0.06
[2,]    0    1 -0.1  0.4    0    0    0    0  0.2   0.0  0.08 -0.02
A matrix matches expected polynomial multiplication result: TRUE 
```

### 3. 完整結果輸出 (Complete Results)

#### F 矩陣 (F Matrix)
```
F Matrix:
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
[1,]    1    0  0.5  0.2    0    0    0    0    0     0     0     0
[2,]    0    1  0.1  0.6    0    0    0    0    0     0     0     0
```

#### A 矩陣 (A Matrix)
```
A Matrix:
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
[1,]    1    0  0.3  0.1    0    0    0    0  0.0   0.2  0.02  0.06
[2,]    0    1 -0.1  0.4    0    0    0    0  0.2   0.0  0.08 -0.02
```

## 驗證結論 (Verification Conclusion)

```
PASS: Complex test matches expected results based on algorithm

Overall validation complete.
```

所有測試均通過，驗證了：
1. R 實現的邏輯與預期的算法相符
2. F 矩陣正確地計算為 [I, FR1]，其中 I 是單位矩陣，FR1 是自回歸係數
3. A 矩陣正確地計算為多項式乘法的結果 (I + AR1*B) * (I + AS1*B^4)
4. 所有矩陣維度正確
5. 結果與理論預期一致

## 技術細節 (Technical Details)

### 算法描述 (Algorithm Description)
`thd2arma` 函數將 THD 模型轉換為簡化形式的 VARMAX 表示法，通過以下步驟：

1. 使用 `thd2arm2` 獲取分解矩陣 (FR, FS, AR, AS, V, G)
2. 對 F(B) = FR(B) * FS(B^s) 進行多項式乘法
3. 對 A(B) = AR(B) * AS(B^s) 進行多項式乘法
4. 構建輸出矩陣 F, A, V, G

### 關鍵參數 (Key Parameters)
- `m`: 內生變量數量 (number of endogenous variables)
- `r`: 外生變量數量 (number of exogenous variables)
- `s`: 季節週期 (seasonal period)
- `p`, `P`: 常規和季節 AR 階數 (regular and seasonal AR orders)
- `q`, `Q`: 常規和季節 MA 階數 (regular and seasonal MA orders)
- `g`: 外生變量滯後階數 (exogenous variable lag order)

### 矩陣維度 (Matrix Dimensions)
- F 和 A 矩陣: m x (m * (k + 1))，其中 k = max(p + P*s, q + Q*s, g)
- G 矩陣: m x (r * (g + 1))
- V 矩陣: m x m (對稱共變異數矩陣 / symmetric covariance matrix)