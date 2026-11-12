"""HTTP下载长度校验；不依赖界面。"""
def check_length(response,received):
    header=getattr(response,'headers',{}).get('Content-Length')
    if header is not None and received!=int(header):raise IOError('下载短读：收到%d字节，HTTP声明%s字节'%(received,header))

