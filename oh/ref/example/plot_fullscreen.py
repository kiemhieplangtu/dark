from matplotlib import pyplot as plt
 

 
# 'Qt4Agg' backend
plt.figure(3)

plt.switch_backend('QT4Agg') #default on my system
plt.plot([1,2,6,4])
mng = plt.get_current_fig_manager()
mng.window.showMaximized()
plt.show()