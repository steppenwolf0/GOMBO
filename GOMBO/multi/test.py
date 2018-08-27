from crossValB import *

test,oneHot_test_labels,valid,oneHot_valid_labels,train,oneHot_train_labels=crossValB_10_Generic(0,10,1)

#8
show(valid[9999])
print(oneHot_valid_labels[9999])
#3
show(train[10])
print(oneHot_train_labels[10])
#0
show(test[750])
print(oneHot_test_labels[750])