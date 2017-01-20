# Create a list.
elements = []

# Append empty lists in first two indexes.
elements.append([])
elements.append([])

# Add elements to empty lists.
elements[0].append([1, 2])
elements[0].append(2)

elements[1].append(3)
print elements[0]
elements[1].append(4)

# Display top-left element.
print(elements[0][0][0])

# Display entire list.
print(elements)


x=0
y=1
z=2
if z>=y and x:
    print "hello"
