def define_network(params):
    import layer_utils
    import metrics
    
    from tensorflow.keras.models import Sequential, Model
    from tensorflow.keras.layers import Dense, Conv1D, Input, BatchNormalization, Activation, concatenate, GlobalAveragePooling1D
#    from tensorflow.keras.utils.np_utils import to_categorical
    from tensorflow.keras.utils import to_categorical
    from tensorflow.keras import optimizers, initializers
    from keras.models import model_from_json
    from keras import layers, models
##    from keras.models import Sequential, Model
##    from keras.layers.core import Dense, Activation
##    from keras.layers import Input, BatchNormalization, Concatenate, Conv1D
##    from tensorflow.keras.utils import to_categorical
##    from keras import optimizers

    # Defining promoter model
    DNA_model = Sequential()
    DNA_model.add(Conv1D(int(params['DNA_n_channel_1st']),
                                int(params['DNA_conv_strides']),
                                input_shape=(None, int(params['n_feature_promoter'])),
                        use_bias=False))
    for k in range(int(params['DNA_n_ConvLayer'])):
        DNA_model.add(Conv1D(int(params['DNA_n_channel_1st']),
                                    int(params['DNA_conv_strides']),
                        use_bias=False))
    DNA_model.add(GlobalAveragePooling1D())

    # Stacking a deep densely-connected network on top
    for k in range(int(params['Last_fullConLayer'])):
        DNA_model.add(Dense(int(params['Last_n_channel'])))
        if params['FullRelu']=='Yes':
            DNA_model.add(Activation('relu'))
    DNA_model.add(Dense(params['n_out'], activation='linear'))

    # Compiling model
    DNA_model.compile(loss='mean_squared_error',
                  optimizer=optimizers.Adam(lr=params['lr']),
                  metrics=[metrics.pcor])

    DNA_model.summary()
    DNA_model.save('test3cs2_model.h5')
    return DNA_model
