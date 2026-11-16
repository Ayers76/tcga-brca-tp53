# 1. 查看当前工作目录
getwd()

# 2. 查看项目文件
list.files()

# 3. 查看 GDCdata 内部结构
list.files("GDCdata", recursive = FALSE)

# 4. 检查两个 RDS
df <- readRDS("df.rds")
results <- readRDS("results.rds")

# 5. 查看对象类型
class(df)
class(results)

# 6. 查看对象结构（限制输出深度）
str(df, max.level = 2)
str(results, max.level = 2)

# 7. 查看文件大小
file.info(c("df.rds", "results.rds"))$size / 1024^2
