
# 【NOTE】
# 参数整理，例如 kpred.mode这种已知所有可能的参数，在函数的参数设置中需要都直接给全。
# 不要加设额外的NULL handle部分了。

# TgView-prep.R DONE

# 382 lines: 这里重新修改一些命名。需要和前面的对齐。考虑修改代码。
# 目标： 避免后续出现意想不到的错误。


# TgView-part2.R [TODO]

# 修改参数plot.x.to.y
# 考虑新增一个参数，让画图的X和Y轴  和 选取需要呈现的 X.to.Y的方向性分离


# TgView-part1.R DONE

# GetResultPieActionMode 需要修改
# 1. color palette的问题，要和action effect统一
# 2. etc

# TgView-part3.R DONE

# GetResultTgSpecificity 需要修改
# 主要是目前的参数太多太乱了，难以理解
# select.genepairs 参数需要重新设计