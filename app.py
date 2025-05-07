# app.py
import os
import pickle
import json
from functools import wraps
from flask import Flask, request, jsonify, render_template, send_from_directory
import lmdb
import csv


app = Flask(__name__)
app.config.from_envvar('FLASK_ENV_FILE', silent=True)

# --------------------------
# 加载基因注释
# --------------------------
GENE_ANNOTATION_PATH = "data/gene_annotation.csv"
gene_id_to_name = {}

with open(GENE_ANNOTATION_PATH) as f:
    reader = csv.reader(f)
    next(reader)  # 跳过标题
    for row in reader:
        gene_id_to_name[row[0]] = row[1]

# --------------------------
# 加载4个 LMDB 环境
# --------------------------
SAMPLES = ["T25LD", "T18DD", "T25DD", "T29DD"]
lmdb_envs = {
    sample: lmdb.open(os.path.join("lmdbs", f"{sample}.lmdb"), readonly=True)
    for sample in SAMPLES
}

# --------------------------
# 访客记录 (SQLite)
# --------------------------
from flask_sqlalchemy import SQLAlchemy
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///visitors.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

class Visitor(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    ip = db.Column(db.String(15), nullable=False)
    country = db.Column(db.String(100))
    visit_time = db.Column(db.DateTime, default=db.func.current_timestamp())

@app.before_request
def track_visitor():
    ip = request.remote_addr
    try:
        res = requests.get(f'https://ipapi.co/{ip}/json/').json()
        country = res.get('country_name')
    except:
        country = None

    visitor = Visitor(ip=ip, country=country)
    db.session.add(visitor)
    db.session.commit()

# --------------------------
# 请求频率限制
# --------------------------
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
# 创建 limiter 实例，指定 key_func 为 IP 地址
limiter = Limiter(
    app=app,
    key_func=get_remote_address,      # 按客户端 IP 做限制
    default_limits=["200 per day", "100 per hour"]  # 可选全局限制
)

# --------------------------
# JSON API 接口
# --------------------------
@app.route('/gene')
@limiter.limit("10 per minute")
def get_gene_expr():
    gene_ids = request.args.get('gene_ids', '').split(',')
    gene_ids = [g.strip() for g in gene_ids if g.strip()]

    if not gene_ids or len(gene_ids) > 100:
        return jsonify({"error": "最多支持100个geneID，用逗号分隔"}), 400

    results = []

    for gene_id in gene_ids:
        if gene_id not in gene_id_to_name:
            continue

        gene_data = {"gene_id": gene_id, "gene_name": gene_id_to_name[gene_id], "expression": {}}

        for sample in SAMPLES:
            with lmdb_envs[sample].begin() as txn:
                data = txn.get(gene_id.encode())
                if data:
                    fpkms = pickle.loads(data)
                    times = sorted(fpkms.keys())
                    gene_data["expression"][sample] = {
                        "time_points": times,
                        "rep1": [fpkms[t]["r1"] for t in times],
                        "rep2": [fpkms[t]["r2"] for t in times],
                        "mean": [fpkms[t]["mean"] for t in times]
                    }

        if gene_data["expression"]:
            results.append(gene_data)

    return jsonify(results)

# --------------------------
# 页面路由
# --------------------------
@app.route('/')
def index():
    return render_template("index.html")

# --------------------------
# 文件下载
# --------------------------
@app.route('/download/<filename>')
def download(filename):
    return send_from_directory("data", filename)

# --------------------------
# 初始化数据库
# --------------------------
with app.app_context():
    db.create_all()

if __name__ == '__main__':
    app.run(debug=False)
