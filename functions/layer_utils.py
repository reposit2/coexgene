#MIT License
#Copyright (c) 2016, Ziga Avsec
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


from tensorflow.keras.layers import AveragePooling2D, AveragePooling1D, Activation
#from keras import backend as K
import tensorflow.keras.backend as K
#import tensorflow.experimental.numpy as tnp
#tnp.experimental_enable_numpy_behavior()
import tensorflow as tf
#tf.compat.v1.disable_v2_behavior()

class GlobalSumPooling1D(AveragePooling1D):
    """Global average pooling operation for temporal data.
    # Note
      - Input shape: 3D tensor with shape: `(batch_size, steps, features)`.
      - Output shape: 2D tensor with shape: `(batch_size, channels)`
    """
    
    def call(self, inputs):
        return K.sum(inputs, axis=1)

class GlobalSumPooling2D(AveragePooling2D):
    """Global sum pooling operation for spatial data.
    # Arguments
        data_format: A string,
            one of `channels_last` (default) or `channels_first`.
            The ordering of the dimensions in the inputs.
            `channels_last` corresponds to inputs with shape
            `(batch, height, width, channels)` while `channels_first`
            corresponds to inputs with shape
            `(batch, channels, height, width)`.
            It defaults to the `image_data_format` value found in your
            Keras config file at `~/.keras/keras.json`.
            If you never set it, then it will be "channels_last".
    # Input shape
        - If `data_format='channels_last'`:
            4D tensor with shape:
            `(batch_size, rows, cols, channels)`
        - If `data_format='channels_first'`:
            4D tensor with shape:
            `(batch_size, channels, rows, cols)`
    # Output shape
        2D tensor with shape:
        `(batch_size, channels)`
    """

    def call(self, inputs):
        if self.data_format == 'channels_last':
            return K.sum(inputs, axis=[1, 2])
        else:
            return K.sum(inputs, axis=[2, 3])

    def get_config(self):
        config = super(ReluFilter, self).get_config()
        config.update({"intercept": self.intercept})
        return config

