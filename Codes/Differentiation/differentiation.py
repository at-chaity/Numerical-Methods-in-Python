
def diff(x,y,p):
    h=x[1]-x[0]
    for j in range(len(x)):
        if x[j]==p:
            i=j
    f_val1=(y[i+1]-y[i])/h
    f_val2=(y[i]-y[i-1])/h
    f_val3=(-3*y[i]+4*y[i+1]-y[i+2])/(2*h)
    f_val4=(y[i+1]-y[i-1])/(2*h)
    f_val5=(y[i-2]-8*y[i-1]+8*y[i+1]-y[i+2])/(12*h)
    return f_val1,f_val2,f_val3,f_val4,f_val5




