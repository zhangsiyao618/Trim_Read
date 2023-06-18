# trim_read README文件

### 说明

trim_read是一个用Python编写的dna测序数据修剪工具，它可以对原始测序数据进行预处理和清理，提高数据的质量和可用性。trim read工具的主要功能是对测序数据进行质量修剪和接头去除，以获得更好的数据解读结果。



### 安装

在运行trim read之前，需要确保Python已经安装在计算机上，并且安装了以下包。

1. `argparse`

   ```bash
   pip install argparse
   ```

2. `matplotlib`

   ```bash
   pip install matplotlib
   ```

3. `numpy`

   ```bash
   pip install numpy
   ```

4. `pandas`

   ```bash
   pip install pandas
   ```



### trim_read的使用方法

进入trim_read文件夹：

```bash
cd 文件夹绝对路径
```

#### 命令行选项

在命令行中运行trim read时，可以使用多个选项，以根据实际情况控制如何进行数据修剪。以下是可选项及其用法：

位置参数：

- 输入`fastq`文件的路径（必选）。

可选参数：

- `-q`或`--fastq`：输出文件格式为fastq（可选）。
- `-a`或`--fasta`：输出文件格式为fasta（可选）。
- `-w`或`--window_size`：设定窗口长度（可选，默认为5）。
- `-t1`或`--threshold1`：设定错误阈值，用以寻找最佳adapter-read匹配时的错误率阈值，同时也是删减长度占比原长度的比例阈值（可选，默认为0.2）。
- `-t2`或`--threshold2`：设定残基质量修剪的质量阈值（可选，默认为10）。
- `l`或`--adapter5_list`：设定`adapter`库文件，文件格式为每行一条序列（可选，默认`adapter`库为adapter5_list.txt）。



#### 示例

下面是针对给定输入数据进行trim read的示例：

输入：

```bash
python trim_read.py .\reads_q64_cut.fastq -a -t1 0.3 -t2 30
```



### 关于反馈和贡献

如果您使用trim read时遇到问题，或者想要为该项目做出改进、提供反馈和建议，请随时与我们联系，我们将不断完善该项目，并得到更好的发展。我们的联系方式：monmon1006@sjtu.edu.cn。