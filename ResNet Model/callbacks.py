import wandb
import tensorflow as tf

# Custom callback that performs batchwise updates

class WandbBatchHistory(tf.keras.callbacks.Callback):
    
    # track the steps and epochs for unique step values
    train_epoch_iter = 0
    train_step_iter = 0
    val_epoch_iter = 0
    val_step_iter = 0
    
    # log training batch metrics
    def on_batch_end(self, batch, logs=None):
        step = batch + self.train_epoch_iter
        wandb.log({'loss':logs['loss']}, step=step, commit=False)
        wandb.log({'accuracy':logs['accuracy']}, step=step)
        self.train_step_iter += 1
        
    # log validation batch metrics
    def on_test_batch_end(self, batch, logs=None):
        step = batch + self.val_epoch_iter
        wandb.log({'val_loss':logs['loss']}, step=step, commit=False)
        wandb.log({'val_accuracy':logs['accuracy']}, step=step)
        self.val_step_iter += 1
    
    # log training epoch final iteration
    def on_epoch_end(self, epoch, logs=None):
        self.train_epoch_iter = self.train_step_iter
    
    # log validation epoch final iteration
    def on_test_epoch_end(self, epoch, logs=None):
        self.val_epoch_iter = self.val_step_iter
        
        