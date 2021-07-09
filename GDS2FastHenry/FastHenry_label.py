# Input file name and desired output name
filename = './example_output.txt'
new_filename = './example_output_labeled.txt'


file = open(filename, 'r')
new_file = open(new_filename, 'w+')

counter = 1
last_line_blank = False
first_coords = ''
for line in file:
    s = str(line)
    # Check for empty lines and repeat the first node of the object if the net line is blank (end of object)
    if s[0] == '\n':
        last_line_blank = True
        new_line = 'N' + str(counter) + first_coords + '\n'
        new_file.write(new_line)
    elif s[0] != '\n':
        s1 = ''
        s2 = ''
        tab = False
        # Loop through each character in line
        for c in s:
            # Store in first string if it's the first number or second string if it's the second number
            if c == '\t':
                tab = True
            elif c != '\t' and tab == False:
                s1 += c
            elif c != '\t' and tab == True:
                s2 += c
        # Create new line with labels
        new_line = 'N' + str(counter) + ' X=' + s1 + '\t' + 'Y=' + s2
        # Save the first node of each object
        if counter == 1 or last_line_blank == True:
            # Change the X coordinate of the first node slightly
            f1 = float(s1) + 0.01
            first_coords = ' X=' + str(f1) + '\t' + 'Y=' + s2
            last_line_blank = False
        new_file.write(new_line)
    counter += 1

# Write the last node of the last object
new_file.write('\nN' + str(counter) + first_coords)
