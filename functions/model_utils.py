def test_prediction2(outloc,
                    best_model,
                    X_promoter_test,
                    Y_test):
    # Inputs
    # outloc: a full path to a directory of the best model
    # best_model: name of the best model
    # X_mRNA_test: mRNA annoation data
    # X_promoter_test: Promoter annoatation data
    # Y_test: Transcriptome data     
    
    # Outputs
    #.{outloc}/test_data/prediction.txt.gz: Predicted gene expression data
    #.{outloc}/test_data/actual.txt.gz: Actual gene expression data
    #.{outloc}/test_data/geneid.txt.gz: Genes in testing data.
    
    import sys
    import os
#    sys.path.append(os.path.abspath("."))
    sys.path.append('./functions/')
    import layer_utils
    import metrics
    from keras.models import load_model
    import data2
#    from . import data
#    import functions.data
    import numpy as np
    import os 
    
    if not os.path.exists(outloc+best_model+'/test_data/'):
        os.makedirs(outloc+"/"+best_model+'/test_data/')
    
    # Load best model
    model = load_model(outloc+best_model+'_model.h5',
                       custom_objects={'pcor': metrics.pcor,
                                      'GlobalSumPooling1D': layer_utils.GlobalSumPooling1D})

    # Batching testing data
    batch_size=128
    test_steps, test_batches = data2.batch_iter2(
                                               X_promoter_test.values[:,1],
                                               Y_test.values[:,1:],
                                               batch_size,
                                               shuffle=False)

    # Making prediction
    pred=[]
    actu=[]
#    print(test_steps)
    for i in range(test_steps):
        a=test_batches.__next__()
        b=model.predict(a[0])
        pred.append(b)
        actu.append(np.vstack(a[1]))

    pred=np.vstack(pred)
    actu=np.vstack(actu)
    
    # Save actual and predicted gene expression
    np.savetxt(outloc+best_model+'/test_data/actual.txt',
               actu, delimiter='\t')
    np.savetxt(outloc+best_model+'/test_data/prediction.txt',
               pred, delimiter='\t')
    X_promoter_test['Name'].to_csv(outloc+best_model+'/test_data/geneid.txt', header=False, index=False, sep='\t')
    
    # gzip text files
    os.system("gzip "+outloc+best_model+'/test_data/actual.txt')
    os.system("gzip "+outloc+best_model+'/test_data/prediction.txt')    
    os.system("gzip "+outloc+best_model+'/test_data/geneid.txt')    
    
def compute_DeepLIFT3(outloc,
             best_model,
             X_promoter_test,
             Y_test):

    # Input
    # outloc: a full path to a directory of the best model
    # best_model: name of the best model
    # X_mRNA_test: mRNA annoation data
    # X_promoter_test: Promoter annoatation data
    # Y_test: Transcriptome data    
    
    # Output
    #.{outloc}/{best_model}/DeepLIFT/DNA_{sample index}.txt.gz: DeepLIFT scores of DNA regulators for each sample. The sample index corresponds to the column index of gene expression data, which starts from 0. DeepLIFT scores of regulators were sparated with commas. Its order is identical to the one appeared in {outloc}/feature_norm_stats.txt
    #.{outloc}/{best_model}/DeepLIFT/RNA_{sample index}}.txt.gz: DeepLIFT scores of RNA regulators for each sample. The sample index corresponds to the column index of gene expression data, which starts from 0. DeepLIFT scores of regulators were sparated with commas. Its order is identical to the one appeared in {outloc}/feature_norm_stats.txt

    import sys
    import os
    sys.path.append('./functions/')
##    import tensorflow as tf
##    tf.compat.v1.disable_v2_behavior()
##    from tensorflow.keras.layers import Dense
##    from tensorflow.keras.models import Model
    import layer_utils
    import metrics
##    from tensorflow.keras.models import load_model
    from keras.models import load_model
    from keras.layers.core import Dense, Activation
    from keras.layers import Input, BatchNormalization, Concatenate, Conv1D
    from keras.models import Model, Sequential
    from keras.models import model_from_json
#    from keras.models import model_from_yaml
    import data
    import data2
    import numpy as np
    import deeplift
    from deeplift.layers import NonlinearMxtsMode
    from deeplift.conversion import kerasapi_conversion as kc
    import tensorflow as tf
    tf.compat.v1.disable_v2_behavior()
    tf.compat.v1.disable_eager_execution()
#    import os

    def my_len(l):
        count = 0
        if isinstance(l, list):
            for v in l:
                count += my_len(v)
            return count
        else:
            return 1

    # Load model
    model = load_model(outloc+best_model+'_model.h5',
                       custom_objects={'pcor': metrics.pcor,
                                       'GlobalSumPooling1D': layer_utils.GlobalSumPooling1D})

    model.summary()

    # Get parameters for the last dens layer
    dens_parameter = model.layers[-1].get_weights()
    
    # Construct a single output model 
    # Adding new layers
    fc = Dense(1,activation='linear',name='out')(model.layers[-2].output)
    new_model = Model(inputs=model.input, outputs=fc)
    new_model = Sequential(layers=new_model.layers)
    new_model2_file = 'tmp_model.h5'

    # Paramter for background distribution
    med_promoter_len=int(np.median([x.shape[1] for x in X_promoter_test.values[:,1]]))
    gene_names_test = X_promoter_test.values[:,0]

    # DeepLIFT score for each sample
    for out_indx in range(dens_parameter[1].shape[0]):
##    for out_indx in range(0,1):
###    for out_indx in range(3,4):
        print(out_indx)

        # Set dens parameters
        new_model.layers[-1].set_weights([dens_parameter[0][:,out_indx:(out_indx+1)],dens_parameter[1][out_indx:(out_indx+1)]])
        new_model.save(new_model2_file)

        if not os.path.exists(outloc+best_model+'/DeepLIFT/'):
            os.makedirs(outloc+"/"+best_model+'/DeepLIFT/')

        # Speficy output file names
        outfile_name_at2=outloc+"/"+best_model+'/DeepLIFT/DNA_'+str(out_indx)+'.txt'

        # Batching testing data
#        batch_size=256*4
##        batch_size=512
        batch_size=2048
#        batch_size=3072
#        batch_size=4096
#        batch_size=8192
        test_steps, test_batches = data2.batch_iter_DeepLIFT2(
                                                            X_promoter_test.values[:,1],
                                                            Y_test.values[:,1:],
                                                            batch_size,
                                                            med_promoter_len,
                                                            shuffle=False)

        for i in range(test_steps):
            
            xs_test,ys_test=next(test_batches)

#            print('tyep(xs_test)',type(xs_test))
#            print('xs_test.shape',xs_test.shape)

            # Reshape background
            xs_background=xs_test * 0.5
###            xs_background=xs_test * 0
#            print('type(xs_background)',type(xs_background))
#            print('xs_background.shape',xs_background.shape)
#            print('out_indx,i',out_indx,i)

            # Compute DeepLIFT scores2
            revealcancel_model = kc.convert_model_from_saved_files(
                                 h5_file=new_model2_file,
                                 nonlinear_mxts_mode=NonlinearMxtsMode.DeepLIFT_GenomicsDefault)
            revealcancel_func = revealcancel_model.get_target_contribs_func(find_scores_layer_idx=0, target_layer_idx=-1)
            shap_values = np.array(revealcancel_func(
                              task_idx=0,
                              input_data_list=[xs_test],
                              input_references_list=[xs_background],
                              batch_size=256,
                              progress_update=None))

            with open(outfile_name_at2, 'a') as f_handle:
                for j in range(shap_values.shape[0]):
                    feature_vector=list(map(str, np.sum(shap_values[j,:,:],axis=0)))
                    out_txt=str(out_indx)+'\t'+gene_names_test[j+i*batch_size]+'\t'+','.join(feature_vector)+'\n'
                    f_handle.write(out_txt)

        # gzip text files
        os.system("gzip "+outfile_name_at2)

