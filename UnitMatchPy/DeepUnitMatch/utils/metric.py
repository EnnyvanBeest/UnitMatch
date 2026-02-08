class AverageMeter(object):
    """Computes and stores the average and current value
       Imported from https://github.com/pytorch/examples/blob/master/imagenet/main.py#L247-L262
    """

    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count

def check_gradients(model, threshold=1e5):
    """
    Checks the gradients of the model's parameters.
    If the norm of the gradients exceeds a certain threshold, training is stopped.

    :param model: The neural network model.
    :param threshold: The threshold for considering a gradient to be exploding.
    :return: A boolean indicating whether training should be stopped.
    """
    for name, param in model.named_parameters():
        if param.grad is not None:
            grad_norm = param.grad.norm().item()
            if grad_norm > threshold:
                print(f"Gradient explosion detected in {name}: {grad_norm}")
                return True  # Indicates training should be stopped
    return False  # Indicates training is fine to continue