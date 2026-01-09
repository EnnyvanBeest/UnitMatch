import logging
import os
import argparse
import tqdm 
import torch
import torch.optim as optim
from torch.utils.data import DataLoader
from torch.utils.tensorboard  import SummaryWriter

from utils import metric
from utils.AE_npdataset import AE_NeuropixelsDataset
from utils.losses import AELoss
from utils.mymodel import *


logger = logging.getLogger(__name__)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def test(model, test_loader, AE_loss, writer):
    model.eval()
    if torch.cuda.is_available():
        model = model.cuda()
    losses = metric.AverageMeter()
    progress_bar = tqdm.tqdm(total=len(test_loader))

    with torch.no_grad():
        for data in test_loader:
            if torch.cuda.is_available():
                data = data.cuda()
            bsz = data.shape[0]

            target = model(data)
            loss = AE_loss(target, data)
            losses.update(loss.item(), bsz)
            progress_bar.update(1)
        progress_bar.close()

    writer.add_scalar('Test/Loss', losses.avg)
    print('Test Loss: %.9f'%(losses.avg))
    return

def validation(epoch, model, val_loader, AE_loss, writer):
    model.eval()
    if torch.cuda.is_available():
        model = model.cuda()

    losses = metric.AverageMeter()
    progress_bar = tqdm.tqdm(total=len(val_loader), desc='Epoch {:3d}'.format(epoch))

    with torch.no_grad():
        for data in val_loader:
            if torch.cuda.is_available():
                data = data.cuda()
            bsz = data.shape[0]
            # Forward pass
            target = model(data)
            loss = AE_loss(target, data)
            losses.update(loss.item(), bsz)
            progress_bar.update(1)
        progress_bar.close()

    writer.add_scalar('Validation/Loss', losses.avg, epoch)
    print('Epoch: %d'%(epoch), 'Val Loss: %.9f'%(losses.avg))
    return

def train(epoch, model, optimizer, train_loader, AE_loss,writer):
    model.train()
    if torch.cuda.is_available():
        model = model.cuda()

    losses = metric.AverageMeter()
    iteration = len(train_loader) * epoch
    
    progress_bar = tqdm.tqdm(total=len(train_loader), desc='Epoch {:3d}'.format(epoch))
    for data in train_loader:
        # data: [bsz,48,2,82], namely [bsz,x_channel, y_channel,time]
        bsz = data.shape[0]
        if torch.cuda.is_available():
            data = data.cuda()
        optimizer.zero_grad()
        target = model(data)
        loss = AE_loss(target, data)
        # update metric
        losses.update(loss.item(), bsz)
        loss.backward()
        # torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=200.0)
        optimizer.step()
        progress_bar.update(1)
        iteration += 1
        if iteration % 50 == 0:
            writer.add_scalar('Train/Loss', losses.avg, iteration)
    
    progress_bar.close()
    print(' Epoch: %d'%(epoch), 'Train Loss: %.9f'%(losses.avg))
    return 

def run(args):
    save_folder = os.path.join('./ModelExp/AE_experiments',  args.exp_name)
    if os.path.exists(save_folder) and not args.cont:
        # Check if the folder already exists and if we are not continuing training
        raise ValueError("Running this experiment will overwrite a previous one. Please choose a different name.")
    ckpt_folder = os.path.join(save_folder,  'ckpt')
    log_folder = os.path.join(save_folder,  'log')
    os.makedirs(ckpt_folder, exist_ok=True)
    os.makedirs(log_folder, exist_ok=True)
    writer = SummaryWriter(log_dir=log_folder)

    np_root = args.train_root
    np_dataset = AE_NeuropixelsDataset(root=np_root, batch_size=args.batchsize)
    print(f"Total size of dataset (train + val + test): {len(np_dataset)}")
    print(f"To open Tensorboard, run this: tensorboard --logdir {os.path.join(os.getcwd(), log_folder)}")
    train_size = int(0.9 * len(np_dataset))
    val_size = int(0.05 * len(np_dataset))
    test_size = len(np_dataset) - train_size - val_size
    train_dataset, val_dataset, test_dataset = torch.utils.data.random_split(np_dataset, [train_size, val_size, test_size])
    train_loader = DataLoader(train_dataset, batch_size=args.batchsize, shuffle=True, num_workers=4)
    val_loader = DataLoader(val_dataset, batch_size=args.batchsize, shuffle=True, num_workers=4)
    test_loader = DataLoader(test_dataset, batch_size=args.batchsize, shuffle=True, num_workers=4)

    model = SpatioTemporalAutoEncoder_V2(n_channel=30,n_time=60,n_output=256).to(device)
    model = model.double()
    encoder = model.encoder

    # AE loss
    AE_Loss = AELoss(lambda1 = 0.,lambda2 = 1.0).to(device)
    optimizer = optim.Adam(model.parameters(), lr=args.lr)

    if args.cont:
        # load latest checkpoint
        ckpt_lst = os.listdir(ckpt_folder)
        ckpt_lst.sort(key=lambda x: int(x.split('_')[-1]))
        read_path = os.path.join(ckpt_folder, ckpt_lst[-1])
        print('load checkpoint from %s'%(read_path))
        checkpoint = torch.load(read_path)
        model.load_state_dict(checkpoint['model'])
        optimizer.load_state_dict(checkpoint['optimizer'])
        encoder.load_state_dict(checkpoint['encoder'])
        start_epoch = checkpoint['epoch'] + 1
    else:
        start_epoch = 0

    if args.total_epoch == 0:
        # don't train, just save the untrained checkpoint
        state = {
            'model': model.state_dict(),
            'encoder': model.encoder.state_dict(),
            'optimizer': optimizer.state_dict(),
            'epoch': 0,
        }
        save_file = os.path.join(ckpt_folder, 'ckpt_epoch_0')
        torch.save(state, save_file)

    for epoch in range(start_epoch, args.total_epoch):
        train(epoch, model, optimizer, train_loader, AE_Loss, writer)
        if epoch % args.save_freq == 0:
            state = {
                'model': model.state_dict(),
                'encoder': model.encoder.state_dict(),
                'optimizer': optimizer.state_dict(),
                'epoch': epoch,
            }
            save_file = os.path.join(ckpt_folder, 'ckpt_epoch_%s'%(str(epoch)))
            torch.save(state, save_file)

        validation(epoch, model, val_loader, AE_Loss, writer)
    
    test(model, test_loader, AE_Loss, writer)
    return 


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--exp_name', '-e', type=str, required=True, help="The checkpoints and logs will be save in ./checkpoint/$EXP_NAME")
    arg_parser.add_argument('--lr', '-l', type=float, default=1e-5, help="Learning rate")
    arg_parser.add_argument('--save_freq', '-s', type=int, default=1, help="frequency of saving model")
    arg_parser.add_argument('--total_epoch', '-t', type=int, default= 300, help="total epoch number for training")
    arg_parser.add_argument('--cont', '-c', action='store_true', help="whether to load saved checkpoints from $EXP_NAME and continue training")
    arg_parser.add_argument('--batchsize', '-b', type=int, default= 32, help="batch size")
    arg_parser.add_argument('--train_root', type=str, default=r"\path\to\your\data", help="root directory of training data")
    args = arg_parser.parse_args()

    run(args)
