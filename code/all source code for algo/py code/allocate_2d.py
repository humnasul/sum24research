def allocate_2d_val(n1, n2):
    p = [None] * n1
    for i in range(n1):
        p[i] = [None] * n2
    return p

'''
arr = allocate_2d_val(5, 4)
arr[0][2] = 6.7
arr[1][3] = 'cat'
print(len(arr[0]))
print(arr[0][2])
print(arr[1][3])
'''