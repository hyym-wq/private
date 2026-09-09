import pandas as pd

file_path = "02_数据/原始数据/附件.xlsx"

df = pd.read_excel(
    file_path,
    sheet_name="男胎检测数据"
)

# -------------------------
# 1. 孕周字符串转成数值
# -------------------------
def week_to_float(value):
    if pd.isna(value):
        return None

    value = str(value).strip().lower()

    if "w+" in value:
        week, day = value.split("w+")
        return float(week) + float(day) / 7

    if "w" in value:
        week = value.replace("w", "")
        return float(week)

    return float(value)
    if pd.isna(value):
        return None

    value = str(value).strip()

    if "w+" in value:
        week, day = value.split("w+")
        return float(week) + float(day) / 7

    if "w" in value:
        week = value.replace("w", "")
        return float(week)

    return float(value)


df["孕周_数值"] = df["检测孕周"].apply(week_to_float)

# -------------------------
# 2. 只提取第一问重要字段
# -------------------------
q1_df = df[
    [
        "孕妇代码",
        "年龄",
        "身高",
        "体重",
        "孕妇BMI",
        "检测孕周",
        "孕周_数值",
        "Y染色体浓度"
    ]
].copy()

# -------------------------
# 3. 查看缺失值
# -------------------------
print("各字段缺失值数量：")
print(q1_df.isna().sum())

print("\n清洗后的前10行：")
print(q1_df.head(10))

# -------------------------
# 4. 保存清洗结果
# -------------------------
output_path = "02_数据/清洗数据/第一问_男胎清洗数据.csv"

q1_df.to_csv(
    output_path,
    index=False,
    encoding="utf-8-sig"
)

print("\n清洗数据已保存到：")
print(output_path)