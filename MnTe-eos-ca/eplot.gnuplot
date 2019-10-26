set terminal postscript
set output "eplot.ps"
#set data style linespoints
set format y "%.4f"
f(x)=a1+a2*x+a3*x**2+a4*x**3+a5*x**4
fit f(x) 'tmp2' via '.fitparam'
plot "tmp2" title "data" w p , f(x) title "polyfit\\_4^{th}order"
#plot "tmp2" title "MnTe-eos-ca"
