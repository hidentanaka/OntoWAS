

import numpy as np
import pandas as pd

data = np.loadtxt('inputfile_OntoT.csv', delimiter=',', dtype=float)
labels = data[:, 0:1] 
features = data[:, 1:] 


Z = features
y_conf = labels.ravel()



from sklearn.decomposition import PCA



pca = PCA(n_components=30)
pca.fit(Z)
X = pca.transform(Z)

import scipy as sp
import pylab as pl

get_ipython().run_line_magic('matplotlib', 'inline')


from limix.ensemble.lmm_forest import  Forest as LMF
from limix.ensemble import lmm_forest_utils as utils

n_samples=965
x = sp.arange(n_samples).reshape(-1,1)




kernel=utils.getQuadraticKernel(x, d=200) + sp.eye(n_samples)*1e-8
n_samples=965


from sklearn.preprocessing import StandardScaler
stdsc = StandardScaler()
(training, test) = utils.crossValidationScheme(2, n_samples)
x_train_std = stdsc.fit_transform(X[training])
x_test_std = stdsc.transform(X[test])

lm_forest = LMF(kernel=kernel[sp.ix_(training, training)])
lm_forest.fit(x_train_std, y_conf[training])
response_tot = lm_forest.predict(x_test_std, kernel[sp.ix_(test,training)])
random_forest = LMF(kernel='iid')
random_forest.fit(x_train_std, y_conf[training])
response_iid = random_forest.predict(x_test_std)

response_fixed = lm_forest.predict(x_test_std)


from matplotlib.backends.backend_pdf import PdfPages 

fig = pl.figure()

pl.plot(x, y_conf, '.7')
pl.plot(x[test], response_tot, 'r-.')
pl.plot(x[test], response_fixed, 'c-.')
pl.plot(x[test], response_iid, 'b-.')
pl.title('prediction')
pl.xlabel('genotype (in decimal encoding)')
pl.ylabel('phenotype')
pl.legend(['fixed effect + confounding', 
           'mixed RF', 'mixed RF (fixed effect)', 'RF'],
           bbox_to_anchor=(1.2, 1.4), ncol=2)
pl.show()

pp = PdfPages('output1_OntT.pdf')

pp.savefig(fig)

pp.close()


from matplotlib.backends.backend_pdf import PdfPages 

fig = pl.figure()


pl.plot(x, y_conf, '.7')
pl.plot(x[test], response_tot, 'r-.')

pl.plot(x[test], response_iid, 'b-.')
pl.title('prediction')
pl.xlabel('genotype (in decimal encoding)')
pl.ylabel('phenotype')
pl.legend(['fixed effect + confounding', 
           'mixed RF', 'RF'],
           bbox_to_anchor=(1.2, 1.4), ncol=2)
pl.show()


pp = PdfPages('output2_OntT.pdf')


pp.savefig(fig)


pp.close()

response_tot_train = lm_forest.predict(X[training], kernel[sp.ix_(training,training)])

response_iid_train = random_forest.predict(X[training])


import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.ensemble import RandomForestRegressor as RFR
import pandas as pd


print('MSE train_mixedRF : %.3f, test_mixedRF : %.3f' % (mean_squared_error(y_conf[training], response_tot_train), mean_squared_error(y_conf[test], response_tot)) )

print('r2 train_mixedRF : %.3f, test_mixedRF : %.3f' % (r2_score(y_conf[training], response_tot_train), r2_score(y_conf[test], response_tot)) )


print('MSE train : %.3f, test : %.3f' % (mean_squared_error(y_conf[training], response_iid_train), mean_squared_error(y_conf[test], response_iid)) )

print('r2 train : %.3f, test : %.3f' % (r2_score(y_conf[training], response_iid_train), r2_score(y_conf[test], response_iid)) )


feature = lm_forest.log_importance
f = pd.DataFrame({'number': range(0, len(feature)),'feature': feature[:]})
f2 = f.sort_values('feature',ascending=False)
pd.set_option('display.max_rows', 1000) 
print(f2)
f2.to_csv("output3.csv")


feature = random_forest.log_importance
f = pd.DataFrame({'number': range(0, len(feature)),'feature': feature[:]})
f2 = f.sort_values('feature',ascending=False)
pd.set_option('display.max_rows', 1000)
print(f2)

f2.to_csv("output4.csv")


import numpy as np
import pandas as pd


data = np.loadtxt('inputfile_gene.csv', delimiter=',', dtype=float, skiprows=1)
labels = data[:, 0:1] 
features = data[:, 1:] 
X = features
y = labels.ravel()


import scipy as sp
import pylab as pl

get_ipython().run_line_magic('matplotlib', 'inline')


from limix.ensemble.lmm_forest import  Forest as LMF
from limix.ensemble import lmm_forest_utils as utils


n_samples=965
x = sp.arange(n_samples).reshape(-1,1)

kernel=utils.getQuadraticKernel(x, d=200) + sp.eye(n_samples)*1e-8

n_samples=965


(training, test) = utils.crossValidationScheme(2, n_samples)

lm_forest = LMF(kernel=kernel[sp.ix_(training, training)])
lm_forest.fit(X[training], y_conf[training])
response_tot = lm_forest.predict(X[test], kernel[sp.ix_(test,training)])

random_forest = LMF(kernel='iid')
random_forest.fit(X[training], y_conf[training])
response_iid = random_forest.predict(X[test])





response_fixed = lm_forest.predict(X[test])





from matplotlib.backends.backend_pdf import PdfPages 

fig = pl.figure()


pl.plot(x, y_conf, '.7')
pl.plot(x[test], response_tot, 'r-.')
pl.plot(x[test], response_fixed, 'c-.')
pl.plot(x[test], response_iid, 'b-.')
pl.title('prediction')
pl.xlabel('genotype (in decimal encoding)')
pl.ylabel('phenotype')
pl.legend(['fixed effect + confounding', 
           'mixed RF', 'mixed RF (fixed effect)', 'RF'],
           bbox_to_anchor=(1.2, 1.4), ncol=2)
pl.show()

# set path
pp = PdfPages('output1_gene.pdf')


pp.savefig(fig)


pp.close()




from matplotlib.backends.backend_pdf import PdfPages 

fig = pl.figure()


pl.plot(x, y_conf, '.7')
pl.plot(x[test], response_tot, 'r-.')
pl.plot(x[test], response_iid, 'b-.')
pl.title('prediction')
pl.xlabel('genotype (in decimal encoding)')
pl.ylabel('phenotype')
pl.legend(['fixed effect + confounding', 
           'mixed RF', 'RF'],
           bbox_to_anchor=(1.2, 1.4), ncol=2)
pl.show()


pp = PdfPages('output2_gene.pdf')


pp.savefig(fig)


pp.close()


response_tot_train = lm_forest.predict(X[training], kernel[sp.ix_(training,training)])

response_iid_train = random_forest.predict(X[training])


import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.ensemble import RandomForestRegressor as RFR
import pandas as pd

print('MSE train_mixedRF : %.3f, test_mixedRF : %.3f' % (mean_squared_error(y_conf[training], response_tot_train), mean_squared_error(y_conf[test], response_tot)) )

print('r2 train_mixedRF : %.3f, test_mixedRF : %.3f' % (r2_score(y_conf[training], response_tot_train), r2_score(y_conf[test], response_tot)) )


print('MSE train : %.3f, test : %.3f' % (mean_squared_error(y_conf[training], response_iid_train), mean_squared_error(y_conf[test], response_iid)) )

print('r2 train : %.3f, test : %.3f' % (r2_score(y_conf[training], response_iid_train), r2_score(y_conf[test], response_iid)) )

feature = lm_forest.log_importance
f = pd.DataFrame({'number': range(0, len(feature)),'feature': feature[:]})
f2 = f.sort_values('feature',ascending=False)
pd.set_option('display.max_rows', 1000) 
print(f2)
f2.to_csv("output3_gene.csv")

feature = random_forest.log_importance
f = pd.DataFrame({'number': range(0, len(feature)),'feature': feature[:]})
f2 = f.sort_values('feature',ascending=False)
pd.set_option('display.max_rows', 1000) 
print(f2)

f2.to_csv("output4_gene.csv")






