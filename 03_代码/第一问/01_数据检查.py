import pandas as pd

# 读取男胎检测数据
file_path = "02_数据/原始数据/附件.xlsx"

df = pd.read_excel(
    file_path,
    sheet_name="男胎检测数据"
)

print("数据行数和列数：")
print(df.shape)

print("\n字段名称：")
print(df.columns.tolist())

print("\n前5行数据：")
print(df.head())

print("\n第一问关键变量：")
print(
    df[
        [
            "孕妇代码",
            "检测孕周",
            "孕妇BMI",
            "Y染色体浓度"
        ]
    ].head(10)
)