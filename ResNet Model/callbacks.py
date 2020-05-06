import wandb
import tensorflow as tf

class WandbBatchHistory(tf.keras.callbacks.Callback,):
        
    def on_train_batch_end(self, batch, logs=None):
        wandb.log({'loss':logs['loss']}, step=batch, commit=False)
        wandb.log({'accuracy':logs['accuracy']}, step=batch)