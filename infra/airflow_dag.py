"""Airflow DAG orchestrating the ETL and analytics pipeline."""
from __future__ import annotations

from datetime import datetime

from airflow import DAG
from airflow.operators.bash import BashOperator


def default_args():
    return {
        "owner": "rx-repurpose",
        "depends_on_past": False,
        "start_date": datetime(2023, 1, 1),
        "retries": 0,
    }


def build_dag() -> DAG:
    with DAG(
        "rx_repurpose_pipeline",
        default_args=default_args(),
        schedule_interval=None,
        catchup=False,
    ) as dag:
        ontology = BashOperator(
            task_id="process_ontology",
            bash_command="python etl/ontology_parser.py",
        )
        tag_ocular = BashOperator(
            task_id="tag_ocular_diseases",
            bash_command="python etl/tag_ocular_diseases.py",
        )
        create_schema = BashOperator(
            task_id="create_graph_schema",
            bash_command="python graph/schema.py --ddl-path graph/schema.sql --manifest-path graph/schema_manifest.json",
        )

        build_graph = BashOperator(
            task_id="build_graph",
            bash_command="python graph/build_graph.py --sqlite-path data/knowledge_graph.sqlite",
        )

        ontology >> tag_ocular >> create_schema >> build_graph
    return dag


dag = build_dag()
