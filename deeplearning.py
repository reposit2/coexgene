# Set input data

# Gene expression matrix (This example uses log2-transformed TPM (Transcripts Per Kilobase Million), but other log2-transformed expression values would work too.) 
deg_data_file = "./data/fourcelltypeswithhffhepg2_fpkm.txt.gz"  # HepG2
#deg_data_file = "./data/fourcelltypeswithhffk562_fpkm.txt.gz"  # K562

# Location of promoter and enhancer features
promoter_data_loc = "./data/"
promoter_annotation_data = ["GTRDE4PROM2006TFM1COR2N890"]   # HepG2
#promoter_annotation_data = ["GTRDE4PROM2006TFM1BCOR2N890"] # K562

# Output directory
outloc='./train_out/3celltypes_hepg2/'  # HepG2
#outloc='./train_out/3celltypes_k562/'  # K562


# Genes used for traning, validation, and testing (A list of gene ids used for training, validation, and testing. These ids should match with gene expression data, RNA features, and promoter fatures)
train_genes = "./data/train.txt.gz"

validate_genes = "./data/validate.txt.gz"

test_genes = "./data/test.txt.gz"

# Location of hyper-parameter
params_loc='./pretrained/Tissue_gene_params7.json' # promoter enhancer output2 one input layer 1311 TFs

import os
#os.environ["CUDA_VISIBLE_DEVICES"]="0"
#import multiprocessing
#process_count = multiprocessing.cpu_count() - 1
import sys
sys.path.append('./functions/')
import data2   # promoter and enhancer
import model_utils
import pandas as pd
import numpy as np
import json
import subprocess
import time
import timeit
import datetime
import tensorflow as tf
from keras.models import model_from_json
#tf.compat.v1.disable_v2_behavior()
#tf.compat.v1.disable_eager_execution()
#os.environ["CUDA_VISIBLE_DEVICES"]="0"
start = time.perf_counter()

def get_h_m_s(td):
    m, s = divmod(td.seconds, 60)
    h, m = divmod(m, 60)
    return h, m, s

def main(params):
    # Descriptions of outputs
    #{time_stamp}_history.json: Model performance at each epoch. Loss and val_loss indicates mean-squred error for traing data and validation data,respectively. Pcor and pcor_loss indicates Pearson's correlation for traing data and validation data,respectively.
    #{time_stamp}_model.h5: The best model among epochs.
    #{time_stamp}_params.json: Hyper parameters.
    #{time_stamp}_test_performance.txt: Model performance of the best model against testing data. The first and the second row indicates mean-squred error and Pearson's correlation, respectively.

#    import keras
    from datetime import datetime
    from tensorflow.python.keras.callbacks import EarlyStopping, ModelCheckpoint
    from tensorflow.python.keras.models import load_model
    import numpy as np
    import metrics
    import layer_utils
    import network
    import tensorflow as tf
#    tf.compat.v1.disable_v2_behavior()
#    tf.compat.v1.disable_eager_execution()

    print(params)
    just_return_model=False

    # model parameters and learning parameters
    max_epoch = 100
    batch_size = 128

    # batch initialization
    train_steps, train_batches = data2.batch_iter2(X_promoter_train.values[:,1],
                                                 Y_train.values[:,1:],
                                                 batch_size,
                                                 shuffle=True)
    valid_steps, valid_batches = data2.batch_iter2(X_promoter_validate.values[:,1],
                                                 Y_validate.values[:,1:],
                                                 batch_size,
                                                 shuffle=True)
    test_steps, test_batches = data2.batch_iter2(X_promoter_test.values[:,1],
                                               Y_test.values[:,1:],
                                               batch_size,
                                               shuffle=True)

    # Paramters for network structure
    params['n_feature_promoter']=X_promoter_train.values[:,1][0].shape[0]
    params['n_out'] = Y_train.values[:,1:].shape[1]

    # Define network structure
    model = network.define_network(params)

    # If you don't need to traning model and just want to have model structure
    if just_return_model:
        return model

    # Set callback functions to early stop training and save the best model so far
    time_stamp=datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    callbacks = [EarlyStopping(monitor='val_loss', patience=10),
                    ModelCheckpoint(outloc+time_stamp+'_model.h5', monitor='val_loss', verbose=0,
                    save_best_only=True,
                    save_weights_only=False,
                    mode='min', period=1)]

    # Optimizing model
    result = model.fit_generator(train_batches, train_steps,
                                 epochs=max_epoch,
                                 validation_data=valid_batches,
                                 validation_steps=valid_steps,
                                 callbacks=callbacks,
                                 max_queue_size=10,
                                 verbose=0)
#                                 use_multiprocessing=True,
#                                 workers=process_count)

    # Test performance
    # Load best model
#    model = load_model(outloc+time_stamp+'_model.h5',
    model = tf.keras.models.load_model(outloc+time_stamp+'_model.h5',
#    model = tf.load_model(outloc+time_stamp+'_model.h5',
##                       custom_objects={'pcor': metrics.pcor})
                       custom_objects={'pcor': metrics.pcor,
                       'GlobalSumPooling1D': layer_utils.GlobalSumPooling1D})
    test_performance= np.array(model.evaluate_generator(test_batches,test_steps))
    np.savetxt(outloc+time_stamp+'_test_performance.txt',
               test_performance,delimiter="\t")

    # Saving optimization history
    with open(outloc+time_stamp+'_history.json', 'w') as f:
        result.history['pcor']=[np.float64(x) for x in result.history['pcor']] # need for json dump
#        result.history['acc']=[np.float64(x) for x in result.history['acc']] # need for json dump
        json.dump(result.history, f)

    # Saving model and learning paramters
    with open(outloc+time_stamp+'_params.json', 'w') as f:
        json.dump(params, f)

    # Return validation loss for model selection
    validation_loss = np.amin(result.history['val_loss'])

#    json_string = model.to_json()
#    open(outloc+time_stamp+'_model.json', 'w').write(json_string)

    return {'loss': validation_loss, 'status': STATUS_OK, 'model': model}


os.makedirs(outloc, exist_ok=True)
shuffle="None"
# Prepare learning data
Y_train, Y_validate, Y_test, X_promoter_train, X_promoter_validate, X_promoter_test = data2.prep_ml_data_split4(
    deg_data_file=deg_data_file,
    promoter_data_loc=promoter_data_loc,
    promoter_annotation_data=promoter_annotation_data,
    train_genes=train_genes,
    validate_genes=validate_genes,
    test_genes=test_genes,
    outloc=outloc,
    shuffle=shuffle)

# Obtain hyper parameters
with open(params_loc) as f:
    params=json.load(f)

# Training model
from hyperopt import STATUS_OK
# train model 10 times
for i in range(10):     # make ten models and select the best model
##for i in range(1):    # make one model
    main(params)

# Descriptions of outputs
#{time_stamp}_history.json: Model performance at each epoch. Loss and val_loss indicates mean-squred error for traing data and validation data,respectively. Pcor and pcor_loss indicates Pearson's correlation for traing data and validation data,respectively.
#{time_stamp}_model.h5: The best model among epochs.
#{time_stamp}_params.json: Hyper parameters.
#{time_stamp}_test_performance.txt: Model performance of the best model against testing data. The first and the second row indicates mean-squred error and Pearson's correlation, respectively.

# Summarizing the result
#! Rscript functions/find_best_model.R "$outloc" &> /dev/null
rsc = f"Rscript functions/find_best_model.R {outloc} &> /dev/null"
print(rsc)
subprocess.call(rsc, shell=True)
while((not os.path.isfile(f"{outloc}summary/loss_parameters.rds")) or (not os.path.isfile(f"{outloc}summary/best_model.txt"))):
    time.sleep(1)

# Descriptions of outputs
#./summary/loss_parameters.rds: Hyper paramters and performance of the 10 training runs.
#./summary/best_model.txt: The best model among the 10 training runs.

# Compute Spearman's correlation between actual and predicted expression for each sample
print('Compute Spearmans correlation between actual and predicted expression for each sample')
# Read the best model
#cwd = os.getcwd()
#rsc = f"cat {cwd}{outloc2}summary/best_model.txt"
rsc = f"cat {outloc}summary/best_model.txt"
print(rsc)
subprocess.call(rsc, shell=True)
while(not os.path.isfile(f"{outloc}summary/best_model.txt")):
    time.sleep(1)
with open(outloc+"summary/best_model.txt") as f:
#with open(cwd+outloc2+"summary/best_model.txt") as f:
    best_model=f.readline().rstrip()
    print(best_model)

# Prediction for test samples with the best model
model_utils.test_prediction2(outloc,
                            best_model,
                            X_promoter_test,
                            Y_test)

# Descriptions of outputs
#.{outloc}/{best_model}/test_data/prediction.txt.gz: Predicted gene expression data
#.{outloc}/{best_model}/test_data/actual.txt.gz: Actual gene expression data
#.{outloc}/{best_model}/test_data/geneid.txt.gz: Genes in testing data.

# Compute prediction accuracy for each sample
rsc = f"Rscript --vanilla --slave functions/calc_performance.R {outloc}{best_model} &> /dev/null"
print(rsc)
subprocess.call(rsc, shell=True)
while(not os.path.isfile(f"{outloc}{best_model}/test_data/cor_tbl.txt")):
    time.sleep(1)

# Descriptions of outputs
#.{outloc}/{best_model}/test_data/cor_tbl.txt: Correlation between acutual and predicted gene expression for each sample
rc = f"{outloc}{best_model}/test_data/cor_tbl.txt"
print(rc)
pd.read_csv(rc,sep="\t")

# Compute average DeepLIFT score for each regulator
print('Compute average DeepLIFT score for each regulator')
# Read the best model
with open(outloc+"summary/best_model.txt") as f:
    best_model=f.readline().rstrip()

# Estimate variable imporance using test samples
model_utils.compute_DeepLIFT3(outloc,
                             best_model,
                             X_promoter_test,
                             Y_test)

# Descriptions of outputs
#.{outloc}/{best_model}/DeepLIFT/DNA_{sample index}.txt.gz: DeepLIFT scores of promoter regulators for each sample. The sample index corresponds to the column index of gene expression data, which starts from 0. DeepLIFT scores of regulators were sparated with commas. Its order is identical to the one appeared in {outloc}/feature_norm_stats.txt

# Concatenate and summarize DeepLIFT scores 
rsc = f"Rscript --vanilla --slave functions/summarize_DeepLIFT2.R {outloc}{best_model} &> /dev/null"
print(rsc)
subprocess.call(rsc, shell=True)
while(not os.path.isfile(f"{outloc}{best_model}/DeepLIFT/promoter_importance_mean.txt")):
    time.sleep(1)
rc = f"{outloc}{best_model}/DeepLIFT/promoter_importance_mean.txt"
print(rc)
pd.read_csv(rc,sep="\t")

process_time = time.perf_counter() - start
td = datetime.timedelta(seconds=process_time)
print('process_time', process_time, 'second')
h, m, s = get_h_m_s(td)
print('process_time', h, 'hour', m, 'minute', s, 'second')

