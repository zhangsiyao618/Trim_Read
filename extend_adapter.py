def reverse_complement_each_line(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    reversed_complement_lines = []
    for line in lines:
        line = line.strip()
        reversed_line = line[::-1]
        base_mapping = str.maketrans('ATCG', 'TAGC')
        complement_line = reversed_line.translate(base_mapping)
        reversed_complement_lines.append(complement_line + '\n')

    # write to file
    with open(file_path, 'a') as f:
        f.write('\n')
        for line in reversed_complement_lines:
            f.write(line)

