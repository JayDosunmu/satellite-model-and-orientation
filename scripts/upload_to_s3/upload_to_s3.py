#! /usr/bin/env python

'''
Given file to upload, S3 target bucket and s3 file name, upload file to the bucket

CLI Args:
  --local_file
  --bucket_name
  --s3_file

ENV variables:
  S3_ACCESS_KEY
  S3_SECRET_KEY
'''

import os
import sys

import argparse
import boto3
from botocore.exceptions import NoCredentialsError

parser = argparse.ArgumentParser('upload_to_s3', description='Script to upload files to AWS S3')
parser.add_argument('--local_file','-l', help='The path of the file to upload', required=True)
parser.add_argument('--bucket_name','-b', help='The target S3 bucket to upload to', required=True)
parser.add_argument('--s3_file','-s', help='The name of the file that will appear in S3',required=True)

ACCESS_KEY = os.environ.get('S3_ACCESS_KEY', '')
SECRET_KEY = os.environ.get('S3_SECRET_KEY', '')

if not ACCESS_KEY:
    print('Unable to read S3 access key from environment. Did you forget to set this?', file=sys.stderr)
    exit(1)

if not SECRET_KEY:
    print('Unable to read S3 secret key from environment. Did you forget to set this?', file=sys.stderr)
    exit(1)

def upload_to_aws(local_file, bucket, s3_file):
    s3 = boto3.client('s3', aws_access_key_id=ACCESS_KEY,
                      aws_secret_access_key=SECRET_KEY)

    try:
        s3.upload_file(local_file, bucket, s3_file)
        print("Upload Successful")
        return True
    except FileNotFoundError:
        print("The file was not found")
        return False
    except NoCredentialsError:
        print("Credentials not available")
        return False

args = parser.parse_args()
local_file = args.local_file
bucket_name = args.bucket_name
s3_file = args.s3_file

uploaded = upload_to_aws(local_file, bucket_name, s3_file)
